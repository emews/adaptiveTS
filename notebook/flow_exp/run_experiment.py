import numpy as np
np.errstate(invalid="ignore",divide="ignore")
import argparse
import pickle
from odin import run_sim

from adaptive_TS import TSBatchBO_CRNGP
import pandas as pd
import os
if os.getcwd().endswith('notebook'):
    os.chdir('flow_exp')

## =================================
## BO experiment settings 

# Create an ArgumentParser object
parser = argparse.ArgumentParser(description = "A script that takes two arguments.")

# Define arguments
parser.add_argument("--exp_path", type = str, 
                    help = "path to the experiment output")
parser.add_argument("--exp_seed", type = int, default = 987654321, 
                    help = "seed to reproduce the whole experiment")
parser.add_argument("--init_npar", type = int, default = 10, 
                    help = "Size of unique parameter settings in initial set of simulations")
parser.add_argument("--nrep", type = int, default = 10, 
                    help = "Number of replicates per unique parameter values to consider")
parser.add_argument("--sim_budget", type = int, default = 700, 
                    help = "Simulation budget")
parser.add_argument("--grid_npar", type = int, default = 100, 
                    help = "Size of unique parameter search grid")
parser.add_argument("--nTS_samp", type = int, default = 30, 
                    help = "New simulations to try at each BO iteration")
# parser$add_argument("--GP_type", type = str, default = "CRNGP", 
#                     help = "GP Surrogate type")
parser.add_argument("--err_sig", type = float, default = 0.5, 
                    help = "gaussian likelihood sd")
parser.add_argument("--prop_sig", type = float, default = 0.3, 
                    help = "proposal distribution sd")
parser.add_argument("--r_path", type = str, 
                    help = "path to r source", default = '.')
parser.add_argument("--exp_id", type = str, 
                    help = "experiment identifier", default = '1')


# Parse the command line arguments
args = parser.parse_args()


exp_path = args.exp_path
exp_seed = args.exp_seed
init_npar = args.init_npar
nrep = args.nrep
sim_budget = args.sim_budget
grid_npar = args.grid_npar
nTS_samp = args.nTS_samp
GP_type = 'crnGP'
err_sig = args.err_sig
prop_sig = args.prop_sig

## Create ground truth
par_true = np.array([0.45, 0.35, 50])
ytrue = pd.read_csv('R_results/ground_truth.csv')
ytrue_full = ytrue['n_SI']
#ytrue = run_sim(par_true)
#ytrue_full = ytrue['n_SI']
#ytrue_vec = ytrue_full
#ytrue_full.to_csv('./y_true_from_py.csv',index=False)

# user input
# exp_path = "/lcrc/project/EMEWS/afadikar/git/adaptiveTS/experiments/exp_1/"
# exp_seed = 23
# init_npar = 5
# nrep = 10
# sim_budget = 200
# grid_npar = 100
# nTS_samp = 30

# # optional input
# err_sig = 0.5
# prop_sig = 0.3

# fixed input for these experiments
sim_func = None
p = 2
ref = -2

## =================================================
## run BO

def load_pickle(path):
    with open(path,'rb') as stream:
        return pickle.load(stream)

if __name__ == "__main__":

    CANDIDATE_SURROGATES = ['crnGP','tabPFN','MixedSingleTaskGP','CustomGP']
    surrogate = ['CustomGP']
    for model in surrogate:
        for i in range(1):
            np.random.seed(exp_seed+i)
            out_adaptive_CRNGP = TSBatchBO_CRNGP(init_npar, nrep, p, sim_budget = sim_budget, 
                                                grid_npar = grid_npar, nTS_samp = nTS_samp, 
                                                ytrue = ytrue_full, ref = ref,
                                                err_sig = err_sig, prop_sig = prop_sig, 
                                                sim_func = sim_func, exp_seed = exp_seed,
                                                adaptive = True,surrogate_func=model)



        out_path = f'py_results/results_{model}_budget_{sim_budget}_rep{i+1}.pkl'

        with open(out_path,'wb') as stream:
            pickle.dump(out_adaptive_CRNGP,stream)