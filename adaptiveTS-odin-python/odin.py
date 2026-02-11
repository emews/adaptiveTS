import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import numba as nb

def SIR(beta = 0.2, gamma = 0.1, S_ini = 1000, I_ini = 10, R_ini = 0,tsteps = 100, seed = None):
    r'''Stochastic SIR simulator based on odin-sim.R'''
    rng = np.random.default_rng(seed)
    S = S_ini
    I = I_ini
    R = R_ini
    T = 0
    n_SI = -1
    n_IR = -1
    out = {
        'S':[S],
        'I':[I],
        'R':[R],
        'n_SI': [n_SI],
        'T': [T]
    }
    N = S + I + R
    for i in range(tsteps):
        S = S - n_SI
        I = I + n_SI - n_IR
        R = R + n_IR
        T = T + 1
        
        p_SI = 1 - np.exp(-beta * I / N)
        p_IR = 1 - np.exp(-gamma)

            

        n_SI = rng.binomial(S,p_SI)
        n_IR = rng.binomial(I,p_IR)
    
        out['S'] += [S]
        out['I'] += [I]
        out['R'] += [R]
        out['T'] += [T]
        out['n_SI'] += [n_SI]

    df = pd.DataFrame(out)
    return df

@nb.njit()
def _SIR_nb(beta = 0.2, gamma = 0.1, S_ini = 10_000, I_ini = 10, R_ini = 0,tsteps = 100, seed = None):
    np.random.seed(seed)
    S = S_ini
    I = I_ini
    R = R_ini
    T = 0
    n_SI = -1
    n_IR = -1
    
    N = S + I + R
    Sout = np.zeros(shape=tsteps, dtype=nb.int64)
    Iout = np.zeros(shape=tsteps, dtype=nb.int64)
    Rout = np.zeros(shape=tsteps, dtype=nb.int64)
    Tout = np.zeros(shape=tsteps, dtype=nb.int64)
    n_SIout = np.zeros(shape=tsteps, dtype=nb.int64)

    for i in range(tsteps):
        S = S - n_SI
        I = I + n_SI - n_IR
        R = R + n_IR
        T = T + 1
        
        p_SI = 1 - np.exp(-beta * I / N)
        p_IR = 1 - np.exp(-gamma)
        n_SI = np.random.binomial(S,p_SI)
        n_IR = np.random.binomial(I,p_IR)

        Sout[i] = S
        Iout[i] = I
        Rout[i] = R
        Tout[i] = T
        n_SIout[i] = n_SI
    return Sout, Iout, Rout, Tout, n_SIout

def SIR_nb(beta = 0.2, gamma = 0.1, S_ini = 10_000, I_ini = 10, R_ini = 0,tsteps = 100, seed = None):
    Sout, Iout, Rout, Tout, n_SIout = _SIR_nb(beta = beta, gamma = gamma,
                                              S_ini = S_ini, I_ini = I_ini, R_ini = R_ini,
                                              tsteps = tsteps, seed = seed)

    out = pd.DataFrame({'S':Sout,
                        'I':Iout,
                        'R':Rout,
                        'T':Tout,
                        'n_SI':n_SIout})
    return out


def run_sim(par,return_err = True, y_true = None):
    beta_range = np.array([0.2, 0.5])
    gamma_range = np.array([0.1, 0.4])
    seed = par[2].astype(int)
    beta = beta_range[0] + par[0] * (beta_range.max() - beta_range.min())
    gamma = gamma_range[0] + par[1]*(gamma_range.max() - gamma_range.min())

    sim = SIR_nb(beta = beta, gamma = gamma, seed = seed,
              S_ini=10_000, I_ini = 10)
    
    if return_err and y_true is not None:
        ans = {'out': ((sim['n_SI']-y_true)**2).sum(),
               'res': sim}
        return ans
    return sim

def run_sim_err(par,y_true):
    ans = run_sim(par,return_err=True,y_true=y_true)
    return ans

def generate_ground_truth():
    par = np.array([0.45,0.35,0])
    sim = run_sim(par,return_err=False)
    sim.to_csv('../ground_truth.csv',index=False)


def test_random_seed():
    seed = 1
    par = np.array([0.45,0.35,0])
    y_true = pd.read_csv('../ground_truth.csv')['n_SI'].values
    out1 = run_sim(par,return_err=True,y_true=y_true)
    out2 = run_sim(par,return_err=True,y_true=y_true)
    assert np.allclose(out1['out'],out2['out'])
    assert np.allclose(out1['res'].values,out2['res'].values)

def test_batch(n=100,plot=0):
    from scipy.stats.qmc import LatinHypercube

    sampler = LatinHypercube(d=2,rng=1)
    X = sampler.random(n=n)
    X = np.hstack([X,np.ones(len(X)).reshape(-1,1)])
    sims = []
    for i, par in enumerate(X):
        sim = run_sim(par)
        sim['sim_id'] = i
        sims.append(sim)
    sims = pd.concat(sims)
    if plot:
        fig, ax = plt.subplots()
        for idx, df in sims.groupby('sim_id'):
            ax.plot(df['T'],df['n_SI'],color='blue',alpha=0.3)
        plt.show()
    return sims

if __name__ == "__main__":
    from time import time
    start = time()
    test_batch(n=10000)
    end = time()
    print(f"{end-start:.02f} seconds")