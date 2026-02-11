import numpy as np
import pandas as pd
import torch
torch.set_default_dtype(torch.float64)
from hetgpy import crnGP, homGP
from scipy.stats import norm
from scipy.stats.qmc import LatinHypercube
from scipy.interpolate import interp1d
from itertools import product
from odin import run_sim, run_sim_err
from flow import train_and_sample
from tabpfn import TabPFNRegressor
from botorch.models import MixedSingleTaskGP
from botorch.models.utils.gpytorch_modules import get_matern_kernel_with_gamma_prior
from botorch.models.gpytorch import GPyTorchModel
from botorch.fit import fit_gpytorch_mll
from gpytorch.mlls import ExactMarginalLogLikelihood
from gpytorch.models import ExactGP
from gpytorch.kernels import IndexKernel
from gpytorch.kernels import MaternKernel
from gpytorch.means import ConstantMean
from gpytorch.distributions import MultivariateNormal
from gpytorch.likelihoods import GaussianLikelihood



class CustomGP(ExactGP,GPyTorchModel):
    def __init__(self, train_inputs, train_targets, likelihood,num_tasks=11):
        super().__init__(train_inputs, train_targets, likelihood)
        self.mean_module = ConstantMean()
        self.task_covariance_module = IndexKernel(num_tasks=num_tasks,rank=10)
        self.covar_module = MaternKernel(nu=2.5,active_dims=[0,1])
    def forward(self, x):
        X_continuous = x[:,0:2]
        X_cat = x[:,-1].int()
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(X_continuous) * self.task_covariance_module(X_cat)
        return MultivariateNormal(mean_x, covar_x)

def tabpfn_onehot(Xsgrid,model):
    df = pd.DataFrame(Xsgrid,columns = [f"X{i}" for i in range(Xsgrid.shape[1])])
    seed_col = df.columns[-1]
    par_cols = df.columns.drop(labels=seed_col)
    seeds = df[seed_col]
    cats = model.seeds
    seeds = pd.Categorical(seeds,categories=cats)

    X_one_hot = np.hstack([df[par_cols].values,
                            pd.get_dummies(seeds).values])
    return X_one_hot


def fit_surrogate(model,X,Y):

    if isinstance(model,crnGP) or model is crnGP:
        model = crnGP()
        model.mle(X,Y,covtype="Matern5_2",known=dict(beta0=0))
        
    elif isinstance(model,TabPFNRegressor) or model is TabPFNRegressor:
        model = TabPFNRegressor(model_path='./tabpfn-v2.5-regressor-v2.5_default.ckpt')
        model.X0 =  X[:,0:X.shape[-1]-1] # works for p = model.X0.shape[1] in crnGP
        model.seeds = sorted(np.unique(X[:,-1].astype(int)))
        
        X_one_hot = tabpfn_onehot(X,model=model)
        model.fit(X_one_hot,Y)
    elif isinstance(model,MixedSingleTaskGP) or model is MixedSingleTaskGP:
        X = torch.from_numpy(X)
        Y = torch.from_numpy(Y).reshape(-1,1)
        model = MixedSingleTaskGP(X, Y, cat_dims=[-1],cont_kernel_factory=get_matern_kernel_with_gamma_prior)
        mll = ExactMarginalLogLikelihood(model.likelihood, model)
        fit_gpytorch_mll(mll)
        model.eval()
        model.X0 =  X[:,0:X.shape[-1]-1] # works for p = model.X0.shape[1] in crnGP
    elif isinstance(model,CustomGP) or model is CustomGP:
        X = torch.from_numpy(X)
        Y = torch.from_numpy(Y) # do not reshape?
        model = CustomGP(train_inputs=X, train_targets=Y,likelihood=GaussianLikelihood())
        mll = ExactMarginalLogLikelihood(likelihood=model.likelihood,model=model)
        fit_gpytorch_mll(mll)
        model.X0 =  X[:,0:X.shape[-1]-1] # works for p = model.X0.shape[1] in crnGP
    else:
        raise ValueError(f"model {model} not recognized")
    return model


def make_grid(*args):
    r'''
    
    Parameters
    ----------
    *args: iterable
        iterable of inputs. must all be ndarray_like and have two dimensions
    
    Returns
    -------
    np.array of cartesian product of all the entries in *args
    '''
    return np.array(list(map(np.concatenate,list(product(*args)))))
def loglik(ysim, ytrue, err_sig):
    r'''
    Likelihood under the normal distribution

    Parameters
    ----------
    ysim: ndarray_like
        simulated outputs
    ytrue: float
        observed output
    err_sig: float
        variance of the normal distribution
    
    Returns
    -------
    log likelihood of the data given the simulation outputs
    '''
    return np.log(norm.pdf(ysim,ytrue,err_sig))

def loglik_design(model,newX,ytrue,err_sig):
    r'''
    Likelihood of the data given the predictions

    Parameters
    ----------
    model: hetgpy model
    newX: ndarray_like
        candidate designs (one row per observation)
    ytrue: float
        ground truth
    err_sig: float
        variance of ytrue distribution
    
    Returns
    -------
    ll: loglikelihoods at newX
    '''
    if isinstance(model,crnGP):
       
        preds = model.predict(x=newX,xprime=newX)
    elif isinstance(model,TabPFNRegressor):
        newX_one_hot = tabpfn_onehot(newX,model=model)
        preds = model.predict(newX_one_hot)
        preds = {'mean':preds}
    elif isinstance(model,MixedSingleTaskGP) or isinstance(model,CustomGP): # treatment is same for botorch models
        with torch.no_grad():
            preds = model.posterior(X = torch.from_numpy(newX),observation_noise=True)
            preds = {'mean': preds.mean.detach().numpy().squeeze()}
    else:
        raise ValueError(f"model {type(model)} not supported")
    ll = loglik(preds['mean'],ytrue,err_sig)
    return ll

def create_grid_flow_CRNGP(nparam, nrep, model, ref, err_sig,lhs_seed = None):
    r'''
    create adaptive grid using normalizing flows

    Parameters
    ----------
    nparam: int
        size of LHS
    nrep: int
        number of seeds
    model: hetgpy model
    ref: int
        ground_truth
    err_sig: float
        variance of ref distribution
    
    Returns
    -------
    new X samples

    
    '''
    p = model.X0.shape[1]

    # lhs design
    s = np.arange(1,nrep+1).reshape(-1,1)
    sampler = LatinHypercube(d=p,rng=np.random.randint(low=1,high=100_000))
    Xgrid_01 = sampler.random(n=nparam)
    Xsgrid_01 = make_grid(Xgrid_01,s)

    # importance weights
    w = loglik_design(newX=Xsgrid_01,model = model, ytrue=ref, err_sig=err_sig)

    w_prime = w - w.max()
    w_norm = np.exp(w_prime) / np.sum(np.exp(w_prime))

    # sample 
    grid_ids = np.random.choice(Xsgrid_01.shape[0],size=Xsgrid_01.shape[0],p=w_norm, replace=True)
    grid_ids = np.unique(grid_ids)
    cols = [f'x{i+1}' for i in range(p)]
    cols.append('r')
    grid_df = pd.DataFrame(Xsgrid_01[grid_ids,:],columns=cols)
    w_filtered = w[grid_ids]

    w_prime = w_filtered - w_filtered.max()
    w_norm  = np.exp(w_prime) / np.sum(np.exp(w_prime))
    grid_df['weight'] = w_norm

    generated_df = train_and_sample(grid_df,nsamples=1000,seeds = s.squeeze())

    return generated_df.values

def TSBatchBO_CRNGP(init_npar,nrep,p,sim_budget,grid_npar,nTS_samp,
                    nTS_iter = None, adaptive = True,ytrue = None,
                    ref = None,err_sig = None, prop_sig = None,
                    sim_func = None, exp_seed = None, surrogate_func = 'crnGP'):
    r'''
    trajectory-oriented optimization via batched Thompson sampling

    Parameters
    ----------
    init_npar: int
        number of initial simulations
    nrep: int
        number of seeds
    p: int
        number of parameters
    sim_budget: int
        number of simulations to run
    grid_npar: int
        size of LHS grid
    nTS_samp: int
        number of Thompson samples
    nTS_iter: None
        not used
    adaptive: bool
        whether to use adaptive grid
    ytrue: float
        ground truth
    ref: float
        reference value for normal distribution
    err_sig: float
        Gaussian likelihood variance
    prop_sig: None
        not used
    sim_func: function
        simulation function
    exp_seed: None
        not used
    Returns
    -------
    dict with entries:
        X_list: simulation parameters run
        y_list: outputs (logged, standardized)
        ynative_list: outputs on native scale
        simout_list: list of full simulation outputs

    '''
    if surrogate_func == 'crnGP':
        surrogate = crnGP
    elif surrogate_func == 'TabPFN':
        surrogate = TabPFNRegressor
    elif surrogate_func == 'MixedSingleTaskGP':
        surrogate = MixedSingleTaskGP
    elif surrogate_func == 'CustomGP':
        surrogate = CustomGP
    else:
        raise ValueError(f"surrogate {surrogate_func} not found")

    if sim_func is None:
        sim_func = run_sim_err
    sampler = LatinHypercube(d=p,rng=np.random.randint(low=1,high=100_000))
    X_01 = sampler.random(n=init_npar)
    s = np.arange(1,nrep+1).reshape(-1,1)
    Xs = make_grid(X_01,s)
    simouts = []
    y = np.full(shape = Xs.shape[0],fill_value=np.nan)
    for i in range(y.shape[0]):
        par = Xs[i,:]
        simout = sim_func(par,y_true=ytrue)
        y[i] = simout['out']
        simouts.append(simout['res'])
    ynative = y.copy()    
    y = np.log(y)
    ym = y.mean(),
    ys = y.std()
    Y = (y - ym) / ys
    
    model = fit_surrogate(surrogate,Xs,Y)

    no_of_sims = len(Y)
    X_list = []
    y_list = []
    simout_list = []
    ynative_list = []

    X_list.append(Xs)
    y_list.append(Y)
    simout_list.append(simouts)
    ynative_list.append(ynative)

    if adaptive:
        tt = 2
        while no_of_sims < sim_budget:
            xnew = next_eval_CRN(model = model,
                                grid_npar=grid_npar,
                                nrep = nrep,
                                ref = ref,
                                err_sig = err_sig,
                                adaptive = True,
                                prop_sig = prop_sig,
                                nTS_samp = nTS_samp,
                                use_flow=True)
            xnew = np.atleast_2d(xnew)
            ynew = np.full(len(xnew),fill_value=np.nan)
            simouts = []
            for i in range(xnew.shape[0]):
                simout = sim_func(xnew[i,:],y_true=ytrue)
                ynew[i] = simout['out']
                simouts.append(simout['res'])
            X_list.append(xnew)
            ynative_list.append(ynew)
            yl = np.log(ynew)
            yl = (yl - ym) / ys
            y_list.append(yl)
            
            simout_list.append(simouts)

            # update surrogate
            Xs = np.concatenate([Xs,xnew],axis=0)
            Y = np.concatenate([Y,yl])
            model = fit_surrogate(surrogate,Xs,Y)
            
            no_of_sims += len(ynew)
            print(f"Iter {tt}, batch size: {len(ynew)}")
            tt += 1
    return {"X_list":X_list, "y_list":y_list, "ynative_list":ynative_list, "simout_list":simout_list}


def next_eval_CRN(model,Xgrid = None,evaluated_ids = None,nTS_samp = None,
                  adaptive = None,grid_npar = None,nrep = None,ref = None,err_sig = None,prop_sig = None ,use_flow = None):
    r'''
    Get next evaluations under the crnGP

    Parameters
    ----------
    model: hetgpy model
    Xgrid: ndarray_like 
        grid points (default None)
    evaluated_ids: ndarray_like
        array of evaluated ids (default None)
    nTS_samp: int
        size of Thompson sample
    adaptive: bool | None
        whether to use adaptive grid
    grid_npar: int
        number of pars in grid
    nrep: int
        number of seeds
    ref: int
        reference value
    err_sig: float
        Gaussian likelihood variance
    prop_sig: float
        not used
    use_flow: bool | None
        whether to use normalizing flow
    
    Returns
    -------
    grid of newX points
    '''
    if adaptive:
        out = adaptive_CRN_TS(model = model,
                            grid_npar = grid_npar,
                            nrep = nrep,
                            ref = ref,
                            err_sig = err_sig,
                            prop_sig = prop_sig,
                            nTS_samp = nTS_samp, 
                            use_flow = use_flow)
    return out
    

def adaptive_CRN_TS(model,grid_npar,nrep,ref,err_sig,prop_sig,nTS_samp,use_flow=True):
    if use_flow:
        Xsgrid = create_grid_flow_CRNGP(nparam=grid_npar,nrep=nrep,
                                      model=model,ref=ref,err_sig=err_sig)
    if isinstance(model,crnGP):
        # predict
        pred = model.predict(x=Xsgrid,xprime=Xsgrid)
        m = pred['mean']
        cov = pred['cov']
        cov = 0.5 * (cov + cov.T)
        eigval, eigvec = np.linalg.eigh(cov)
        if (eigval < 0).any():
            # threshold eigvalues to zero and re-calculate cov
            eigval[eigval < 0] = 0
            cov = eigvec @ np.diag(eigval) @ eigvec.T
        tTS = np.random.multivariate_normal(m,cov, size = nTS_samp)
        best_ids = np.argmin(tTS,axis=1)
        best_ids = np.unique(best_ids,axis=0)
    elif isinstance(model,MixedSingleTaskGP) or isinstance(model,CustomGP):
        with torch.no_grad():
            preds = model.posterior(X= torch.from_numpy(Xsgrid),observation_noise=True)
            tTS = preds.sample(sample_shape=torch.Size([nTS_samp]))
            best_ids = np.argmin(tTS,axis=1).squeeze()
            best_ids = np.unique(best_ids,axis=0)

    elif isinstance(model,TabPFNRegressor):
        qs = np.linspace(0.01,0.99,num=30)
        
        X_one_hot = tabpfn_onehot(Xsgrid=Xsgrid,model=model)
        pred = model.predict(X=X_one_hot,quantiles=qs,output_type="full")
        quantiles = np.array(pred['quantiles'])
        idxs = []
        for i in range(nTS_samp):
            # generate a sample from the predictive dist
            out = []
            for x in quantiles.T:
                inverse_cdf = interp1d(qs, x, kind='linear', bounds_error=False, fill_value=(x[0], x[-1]))
                uniform = np.random.uniform(size = 1)
                sample = inverse_cdf(uniform)
                out.append(sample)
            out = np.array(out)
            idx = out.argmin()
            idxs.append(idx)
        best_ids = np.unique(idxs)
    else:
        raise ValueError(f"model {type(model)} not supported")
            

    return Xsgrid[best_ids,:]
    