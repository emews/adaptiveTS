library(hetGP, lib.loc = "/lcrc/project/EMEWS/bebop/R-packages/")
library(scales)
library(data.table)
library(lhs)
library(ggplot2)
library(mvtnorm)
library(argparse)

source("odin-sim.R")
source("plots.R")


## Create ground truth

par_true <- c(0.45, 0.35, 50)
xtrue <- par_true[1:2]
strue <- par_true[3]
ytrue <- run_sim(c(xtrue, strue))
ytrue_full <- ytrue$n_SI
ytrue_vec <- ytrue_full


## =================================
## BO experiment settings 

# Create an ArgumentParser object
parser <- ArgumentParser(description = "A script that takes two arguments.")

# Define arguments
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
parser$add_argument("--GP_type", type = "character", default = "CRNGP", 
                    help = "GP Surrogate type")
parser$add_argument("--err_sig", type = "numeric", default = 0.5, 
                    help = "gaussian likelihood sd")
parser$add_argument("--prop_sig", type = "numeric", default = 0.3, 
                    help = "proposal distribution sd")


# Parse the command line arguments
args <- parser$parse_args()

# exp_seed <- args$exp_seed
# init_npar <- args$init_npar
# nrep <- args$nrep
# sim_budget <- args$sim_budget
# grid_npar <- args$grid_npar
# nTS_samp <- args$nTS_samp
# GP_type <- args$GP_type
# err_sig <- args$err_sig
# prop_sig <- args$prop_sig


# user input
exp_seed <- 23
init_npar <- 5
nrep <- 10
sim_budget <- 200
grid_npar <- 100
nTS_samp <- 30
GP_type <- "CRNGP"

# optional input
err_sig <- 0.5
prop_sig <- 0.3

# fixed input for these experiments
sim_func <- run_sim_err
p <- 2
ref <- -2

## =================================================
## run BO

out_fixed <- TSBatchBO(init_npar, nrep, p, sim_budget = sim_budget, 
                       grid_npar = grid_npar, nTS_samp = nTS_samp, 
                       ytrue = ytrue_full, GP_type = GP_type, ref = ref,
                       err_sig = err_sig, prop_sig = prop_sig, 
                       sim_func = run_sim_err, exp_seed = exp_seed, covtype = "Matern5_2",
                       adaptive = F)

out_adaptive <- TSBatchBO(init_npar, nrep, p, sim_budget = sim_budget, 
                          grid_npar = grid_npar, nTS_samp = nTS_samp, 
                          ytrue = ytrue_full, GP_type = GP_type, ref = ref,
                          err_sig = err_sig, prop_sig = prop_sig, 
                          sim_func = run_sim_err, exp_seed = exp_seed, covtype = "Matern5_2",
                          adaptive = T)

## =================================================
## collect and format results

x_fixed_mat <- do.call(rbind, out_fixed[[1]])
x_adaptive_mat <- do.call(rbind, out_adaptive[[1]])

y_fixed_vec <- unlist(out_fixed[[2]])
y_adaptive_vec <- unlist(out_adaptive[[2]])

## run simulations

y_fixed_traj <- matrix(NA, nrow = 101, ncol = nrow(x_fixed_mat))
for (ii in 1:nrow(x_fixed_mat)){
  o <- run_sim(x_fixed_mat[ii, ])
  y_fixed_traj[, ii] <- o$n_SI
}

y_adaptive_traj <- matrix(NA, nrow = 101, ncol = nrow(x_adaptive_mat))
for (ii in 1:nrow(x_fixed_mat)){
  o <- run_sim(x_adaptive_mat[ii, ])
  y_adaptive_traj[, ii] <- o$n_SI
}

# run simulations and calculate sum squared errors

y_fixed_native <- y_adaptive_native <- rep(NA, sim_budget)
for(ii in 1:sim_budget){
  y_fixed_native[ii] <-run_sim_err(x_fixed_mat[ii, ])
  y_adaptive_native[ii] <- run_sim_err(x_adaptive_mat[ii, ])
}


## Dense evaluations for comparison
nparam_dense <- 50
X_01 <- as.matrix(expand.grid(seq(0, 1, length.out = nparam_dense), 
                              seq(0, 1, length.out = nparam_dense)))
s <- 1:nrep
Xs_01 <- cbind(X_01[rep(1:nparam_dense^2, each = nrep), ], 
               rep(s, nparam_dense^2))

y_dense <- rep(NA, nparam_dense^2*nrep)
k <- 1
for (ii in 1:nparam_dense^2){
  for (seed in s){
    y_dense[k] <- run_sim_err(Xs_01[k, ])
    k <- k + 1
  }
}
d_dense <- data.frame(cbind(Xs_01, y_dense))
colnames(d_dense) <- c("beta", "gamma", "seed", "y")

## ================================================
## plot

p_pairs <- pairs_plot(x_fixed_mat, x_adaptive_mat)
p_pairs

for(ii in c(10, 30, 50, 80, 100)){
  print(plot_nbest_traj(nbest = ii, x_fixed_mat, x_adaptive_mat, y_fixed_vec, y_adaptive_vec, yobs = ytrue_vec))
  print(err_hist(nbest = ii, y_fixed_native, y_adaptive_native))
}


nbests <- seq(20, sim_budget, length.out = 10)
par(mfrow = c(5, 2))
for(nbest in nbests){
  plot_pct_traj(nbest, y_fixed_native, y_adaptive_native)
}

contour_compare(nbest = 50, 1, x_fixed_mat, x_adaptive_mat, y_fixed_vec, y_adaptive_vec, 
                d_dense[d_dense$seed == 1, ])
