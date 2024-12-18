library(hetGP, lib.loc = "/lcrc/project/EMEWS/bebop/R-packages/")
library(scales)
library(data.table)
library(lhs)
library(ggplot2)
library(mvtnorm)
library(argparse)


## =================================
## BO experiment settings 

# Create an ArgumentParser object
parser <- ArgumentParser(description = "A script that takes two arguments.")

# Define arguments
parser$add_argument("--exp_path", type = "character", 
                    help = "path to the experiment output")
parser$add_argument("--exp_seed", type = "integer", default = 987654321, 
                    help = "seed to reproduce the whole experiment")
parser$add_argument("--init_npar", type = "integer", default = 10, 
                    help = "Size of unique parameter settings in initial set of simulations")
parser$add_argument("--nrep", type = "integer", default = 10, 
                    help = "Number of replicates per unique parameter values to consider")
parser$add_argument("--sim_budget", type = "integer", default = 200, 
                    help = "Simulation budget")
parser$add_argument("--grid_npar", type = "integer", default = 100, 
                    help = "Size of unique parameter search grid")
parser$add_argument("--nTS_samp", type = "integer", default = 30, 
                    help = "New simulations to try at each BO iteration")
# parser$add_argument("--GP_type", type = "character", default = "CRNGP", 
#                     help = "GP Surrogate type")
parser$add_argument("--err_sig", type = "numeric", default = 0.5, 
                    help = "gaussian likelihood sd")
parser$add_argument("--prop_sig", type = "numeric", default = 0.3, 
                    help = "proposal distribution sd")
parser$add_argument("--r_path", type = "character", 
                    help = "path to r source", default = '.')
parser$add_argument("--exp_id", type = "character", 
                    help = "experiment identifier", default = '1')


# Parse the command line arguments
args <- parser$parse_args()

source(paste0(args$r_path, "/odin-sim.R"))
source(paste0(args$r_path, "/adaptive_TS.R"))
source(paste0(args$r_path, "/adaptive_TS_nonCRN.R"))
source(paste0(args$r_path, "/plots.R"))

exp_path <- args$exp_path
exp_seed <- args$exp_seed
init_npar <- args$init_npar
nrep <- args$nrep
sim_budget <- args$sim_budget
grid_npar <- args$grid_npar
nTS_samp <- args$nTS_samp
GP_type <- args$GP_type
err_sig <- args$err_sig
prop_sig <- args$prop_sig

## Create ground truth
par_true <- c(0.45, 0.35, 50)
xtrue <- par_true[1:2]
strue <- par_true[3]
ytrue <- run_sim(c(xtrue, strue))
ytrue_full <- ytrue$n_SI
ytrue_vec <- ytrue_full


# user input
# exp_path <- "/lcrc/project/EMEWS/afadikar/git/adaptiveTS/experiments/exp_1/"
# exp_seed <- 23
# init_npar <- 5
# nrep <- 10
# sim_budget <- 200
# grid_npar <- 100
# nTS_samp <- 30

# # optional input
# err_sig <- 0.5
# prop_sig <- 0.3

# fixed input for these experiments
sim_func <- run_sim_err
p <- 2
ref <- -2

## =================================================
## run BO

out_fixed_CRNGP <- TSBatchBO_CRNGP(init_npar, nrep, p, sim_budget = sim_budget, 
                                   grid_npar = grid_npar, nTS_samp = nTS_samp, 
                                   ytrue = ytrue_full, ref = ref,
                                   err_sig = err_sig, prop_sig = prop_sig, 
                                   sim_func = run_sim_err, exp_seed = exp_seed, covtype = "Matern5_2",
                                   adaptive = F)

out_adaptive_CRNGP <- TSBatchBO_CRNGP(init_npar, nrep, p, sim_budget = sim_budget, 
                                      grid_npar = grid_npar, nTS_samp = nTS_samp, 
                                      ytrue = ytrue_full, ref = ref,
                                      err_sig = err_sig, prop_sig = prop_sig, 
                                      sim_func = run_sim_err, exp_seed = exp_seed, covtype = "Matern5_2",
                                      adaptive = T)

out_fixed_hetGP <- TSBatchBO_hetGP(init_npar, nrep, p, sim_budget = sim_budget, 
                                   grid_npar = grid_npar, nTS_samp = nTS_samp, 
                                   ytrue = ytrue_full, ref = ref,
                                   err_sig = err_sig, prop_sig = prop_sig, 
                                   sim_func = run_sim_err, exp_seed = exp_seed, covtype = "Matern5_2",
                                   adaptive = F)

out_adaptive_hetGP <- TSBatchBO_hetGP(init_npar, nrep, p, sim_budget = sim_budget, 
                                      grid_npar = grid_npar, nTS_samp = nTS_samp, 
                                      ytrue = ytrue_full, ref = ref,
                                      err_sig = err_sig, prop_sig = prop_sig, 
                                      sim_func = run_sim_err, exp_seed = exp_seed, covtype = "Matern5_2",
                                      adaptive = T)

save.image(file = paste0(exp_path, "/", args$exp_id, "_", "out.RData"))