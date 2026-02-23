'''
emulator.py

This module stores wrapper code for various emulators
'''
from copy import deepcopy
import torch
import numpy as np
import pandas as pd
# hetgpy inputs
from hetgpy import crnGP, homGP, hetGP
# botorch
from botorch.models.gpytorch import GPyTorchModel
from botorch.fit import fit_gpytorch_mll
from gpytorch.mlls import ExactMarginalLogLikelihood
# gpytorch
from gpytorch.models import ExactGP
from gpytorch.kernels import IndexKernel
from gpytorch.kernels import MaternKernel
from gpytorch.kernels import RBFKernel
from gpytorch.means import ConstantMean
from gpytorch.distributions import MultivariateNormal
from gpytorch.likelihoods import GaussianLikelihood
# tabPFN
from tabpfn import TabPFNRegressor

__all__ = ['BaseEmulator','CategoricalSeedGPEmulator','CategoricalSeedGPEmulatorOnesDiagonal'
           'crnGPEmulator','hetGPEmulator','homGPEmulator','tabPFNEmulator']

def _check_seed(X):
    '''
    Check that last column of X is integer
    '''
    if isinstance(X,torch.Tensor):
        if not (X[:,-1] == X[:,-1].int()).all(): 
            raise RuntimeError("Last column of X must specify integer random seed")
    elif isinstance(X,np.ndarray):
        if not (X[:,-1] == X[:,-1].astype(int)).all():
                raise RuntimeError("Last column of X must specify integer random seed")
    return

class BaseEmulator:
    '''
    Base class to create an emulator. Requires implemening a fit, predict, and sample method
    '''
    def __init__(self):
        return

    def fit(self,X,Y):
        raise NotImplementedError()
    def predict(self):
        raise NotImplementedError()
    def sample(self):
        raise NotImplementedError()
    def copy(self):
        '''
        Return a copy of the model
        '''
        newmodel = deepcopy(self)
        return newmodel

class _CategoricalSeededGP(ExactGP,GPyTorchModel):
    def __init__(self, train_inputs, train_targets, likelihood = GaussianLikelihood(),covtype="Matern5_2",num_tasks=10,rank = 10):
        super().__init__(train_inputs, train_targets, likelihood)
        _check_seed(train_inputs)
        self.mean_module = ConstantMean()
        self.task_covariance_module = IndexKernel(num_tasks=num_tasks,rank=rank)
        self.active_dims = np.arange(train_inputs.shape[1] - 1).tolist()
        ard_num_dims = len(self.active_dims)
        # set covariance type
        if covtype == 'Matern5_2':
            covar = MaternKernel(nu=2.5,active_dims=self.active_dims,ard_num_dims=ard_num_dims)
        elif covtype == 'Gaussian':
            covar = RBFKernel(active_dims=self.active_dims,ard_num_dims=ard_num_dims)
        self.covar_module = covar
    def forward(self, x):
        X_continuous = x[:,self.active_dims]
        X_cat = (x[:,-1] - 1).int() # categoricals must map to {0, ..., (N - 1)}
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(X_continuous) * self.task_covariance_module(X_cat)
        return MultivariateNormal(mean_x, covar_x)
    
class CategoricalSeedGPEmulator(BaseEmulator):
    def __init__(self, likelihood = GaussianLikelihood(),covtype="Matern5_2",num_tasks = None,rank = None):
        super().__init__()
        self.likelihood = likelihood
        self.num_tasks = num_tasks
        self.rank = rank
        self.covtype = covtype
        self.model = None
    def fit(self,X,Y):
        if self.model is not None:
            del self.model
            self.model = None
        if isinstance(X,np.ndarray):
            X = torch.from_numpy(X)
        if isinstance(Y,np.ndarray):
            Y = torch.from_numpy(Y)

        model = _CategoricalSeededGP(train_inputs=X, train_targets=Y,likelihood=self.likelihood,
                                     covtype=self.covtype,
                                     num_tasks=self.num_tasks, rank = self.rank)
        model.train()
        mll = ExactMarginalLogLikelihood(likelihood=model.likelihood,model=model)
        fit_gpytorch_mll(mll)
        
        model.eval()
        self.model = model

    def predict(self,X):
        with torch.no_grad():
            if isinstance(X,np.ndarray):
                X = torch.from_numpy(X)
            post = self.model.posterior(X, observation_noise=True)
            preds = {'mean':post.mean.detach().numpy().squeeze()}
        return preds

    def sample(self,X,size=1):
        with torch.no_grad():
            if isinstance(X,np.ndarray):
                X = torch.from_numpy(X)
            preds = self.model.posterior(X=X,observation_noise=True)

        return preds.sample(sample_shape=torch.Size([size])).squeeze()

class _CustomIndexKernel(IndexKernel):
    def __init__(self, num_tasks, rank = 1, prior = None, var_constraint = None, **kwargs):
        super().__init__(num_tasks, rank, prior, var_constraint, **kwargs)

    def _eval_covar_matrix(self):
        cf = self.covar_factor
        BBT = cf @ cf.transpose(-1, -2)
        diag = BBT.diag()
        out = BBT / torch.sqrt(diag.unsqueeze(-1) * diag.unsqueeze(-2))
        return out + torch.diag_embed(self.var)
        

class _CategoricalSeededGPOnesDiagonal(ExactGP,GPyTorchModel):
    def __init__(self, train_inputs, train_targets, likelihood,covtype="Matern5_2", num_tasks = None,rank = None):
        super().__init__(train_inputs, train_targets, likelihood)
        _check_seed(train_inputs)
        self.mean_module = ConstantMean()
        self.task_covariance_module = _CustomIndexKernel(num_tasks=num_tasks,rank=rank)
        self.active_dims = np.arange(train_inputs.shape[1] - 1).tolist()
        ard_num_dims = len(self.active_dims)
        # set covariance type
        if covtype == 'Matern5_2':
            covar = MaternKernel(nu=2.5,active_dims=self.active_dims,ard_num_dims=ard_num_dims)
        elif covtype == 'Gaussian':
            covar = RBFKernel(active_dims=self.active_dims,ard_num_dims=ard_num_dims)
        self.covar_module = covar

    def forward(self, x):
        X_continuous = x[:,self.active_dims]
        X_cat = (x[:,-1] - 1).int() # categoricals must map to {0, ..., (N - 1)}
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(X_continuous) * self.task_covariance_module(X_cat)
        return MultivariateNormal(mean_x, covar_x)

class CategoricalSeedGPEmulatorOnesDiagonal(BaseEmulator):
    def __init__(self,likelihood = GaussianLikelihood(),covtype = "Matern5_2",num_tasks = None, rank = None):
        super().__init__()
        self.likelihood = likelihood
        self.num_tasks = num_tasks
        self.rank = rank
        self.covtype = covtype
    def fit(self,X,Y):
        self.model = None
        if isinstance(X,np.ndarray):
            X = torch.from_numpy(X)
        if isinstance(Y,np.ndarray):
            Y = torch.from_numpy(Y)

        model = _CategoricalSeededGPOnesDiagonal(train_inputs = X, train_targets = Y,
                                                 likelihood = self.likelihood,
                                                 covtype = self.covtype,
                                                 num_tasks = self.num_tasks, rank = self.rank)
        model.train()
        mll = ExactMarginalLogLikelihood(likelihood=model.likelihood,model=model)
        fit_gpytorch_mll(mll)
        
        model.eval()
        self.model = model

    def predict(self,X):
        with torch.no_grad():
            if isinstance(X,np.ndarray):
                X = torch.from_numpy(X)
            post = self.model.posterior(X, observation_noise=True)
            preds = {'mean':post.mean.detach().numpy().squeeze()}
        return preds

    def sample(self,X,size=1):
        with torch.no_grad():
            if isinstance(X,np.ndarray):
                X = torch.from_numpy(X)
            preds = self.model.posterior(X=X,observation_noise=True)

        return preds.sample(sample_shape=torch.Size([size])).squeeze()






class _hetGPyEmulator(BaseEmulator):
    def __init__(self,covtype = "Gaussian", lower = None, upper = None, known = {}, init = {}):
        self.covtype = covtype
        self.lower = lower
        self.upper = upper
        self.known = known
        self.init = init

        self.base_model = None

    def fit(self,X,Y):
        self.model = self.base_model() # instantiate model object
        self.model.mle(X=X,Z=Y,covtype = self.covtype, lower = self.lower, upper = self.upper)
    
    def predict(self,X,xprime = None):
        preds = self.model.predict(x=X,xprime=xprime)
        return preds
    def sample(self,X,size = 1):
        preds = self.model.predict(x=X,xprime=X)
        # check covariance matrix
        m = preds['mean']
        cov = preds['cov']
        cov = 0.5 * (cov + cov.T)
        eigval, eigvec = np.linalg.eigh(cov)
        if (eigval < 0).any():
            # threshold eigvalues to zero and re-calculate covariance matrix
            eigval[eigval < 0] = 0
            cov = eigvec @ np.diag(eigval) @ eigvec.T
        return np.random.multivariate_normal(m,cov, size = size)
    
class crnGPEmulator(_hetGPyEmulator):
    def __init__(self,covtype="Gaussian", lower=None, upper=None, known={}, init={}):
        super().__init__(covtype, lower, upper, known, init)
        self.base_model = crnGP

class hetGPEmulator(_hetGPyEmulator):
    def __init__(self,covtype="Gaussian", lower=None, upper=None, known={}, init={}):
        super().__init__(covtype, lower, upper, known, init)
        self.base_model = hetGP

class homGPEmulator(_hetGPyEmulator):
    def __init__(self,covtype="Gaussian", lower=None, upper=None, known={}, init={}):
        super().__init__(covtype, lower, upper, known, init)
        self.base_model = homGP

class tabPFNEmulator(BaseEmulator):
    '''
    Wrapper for the tabPFN model with a categorical input. Categorical input must be the last column?
    '''
    def __init__(self,model_path,seed_col = -1):
        super().__init__()
        self.seed_col = seed_col
        self.model = TabPFNRegressor(model_path=model_path,inference_precision=torch.float32,device='cpu')
    def onehot(self,Xsgrid):
        df = pd.DataFrame(Xsgrid,columns = [f"X{i}" for i in range(Xsgrid.shape[1])])
        seed_col = df.columns[-1]
        par_cols = df.columns.drop(labels=seed_col)
        seeds = df[seed_col]
        cats = self.seeds
        seeds = pd.Categorical(seeds,categories=cats)

        X_one_hot = np.hstack([df[par_cols].values,
                                pd.get_dummies(seeds).values])
        return X_one_hot
    def fit(self,X,Y):
        _check_seed(X)
        self.seeds = np.sort(np.unique(X[:,self.seed_col].astype(int)))
        X_onehot = self.onehot(Xsgrid=X)
        X_onehot = torch.from_numpy(X_onehot.astype(np.float32))
        Y = torch.from_numpy(Y.astype(np.float32))
        self.model.fit(X_onehot,Y)
        return
    def predict(self,X,qs = None):
        if qs is None:
            qs = np.linspace(0.01,0.99,num=30)
        self.qs = qs
        X_one_hot = self.onehot(Xsgrid=X)
        preds = self.model.predict(X=X_one_hot,quantiles=qs,output_type="full")
        
        return preds
    def sample(self,X,size=1,qs = None, t=1.0):
        preds = self.predict(X=X,qs=qs)
        criterion, logits = preds['criterion'], preds['logits']
        samples = np.concat([
            criterion.sample(logits,t=t).detach().numpy().reshape(1,-1) 
                for i in range(size)
            ],axis=0)
        return samples


if __name__ == "__main__":
    from hetgpy.example_data import mcycle
    m = mcycle()
    X, Y = m['times'], m['accel']
    Xgrid = np.linspace(X.min(),X.max(),num=200).reshape(-1,1)
    Xgrid = np.concat([Xgrid,
                      np.random.default_rng(1).integers(low=1,high=3,size=len(Xgrid)).reshape(-1,1)],axis=1)
    tab = tabPFNEmulator(model_path='tabpfn-v2.5-regressor-v2.5_default.ckpt')
    tab.fit()
