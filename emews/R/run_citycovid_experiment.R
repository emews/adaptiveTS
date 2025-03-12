library(hetGP, lib.loc = "/lcrc/project/EMEWS/bebop/R-packages/")
library(scales)
library(data.table)
library(lhs)
library(ggplot2)
library(mvtnorm)
library(argparse)
library(reticulate)
library(EQ.SQL)

source(paste0(args$r_path, "/adaptive_TS_citycovid.R"))
gt_h_file <- "../experiments/data/dt.chicago.hosp.csv"

## =================================
## BO experiment settings 


# hard-coding these for now, can make input if necessary
exp_design <- list(
  sim_budget = 2000, # Total simulation budget
  init_npar = 30,    # Number of initial design points
  nrep = 30,         # Number of replications
  grid_npar = 10,    # Number of candidate points for grid
  nTS_samp = 100,    # Number of Thompson samples
  p = 3             # Dimension of the input space
)

exp_seed <- 1 # Is this correct to set here?

## =================================
## EMEWS setup

config_file =  "algo_cfg.yaml"
params = parse_yaml_cfg(config_file)
eqsql <- init_eqsql(python_path = params$python_path)

eqsql$db_tools$stop_db(params$db_path)

db_started <- FALSE
pool <- NULL
task_queue <- NULL

eqsql <- init_eqsql(python_path = params$python_path)

eqsql$db_tools$start_db(params$db_path)
db_started <- TRUE

task_queue <- init_task_queue(eqsql, params$db_host, params$db_user, params$db_port,
                              params$db_name)

if (!task_queue$are_queues_empty()) {
  print("WARNING: task input / output queues are not empty. Aborting run")
  task_queue$clear_queues()
}
pool_params <- eqsql$worker_pool$cfg_file_to_dict(params$pool_cfg_file)
params$worker_pool_id <- 'adaptive_ts'
exp_id <- 1
pool <- eqsql$worker_pool$start_local_pool(params$worker_pool_id, params$pool_launch_script, exp_id, pool_params)
task_type <- params$task_type

## =================================
## run adaptive TS

out <- runAdaptiveTS(exp_design, task_queue=task_queue, exp_id=exp_id, task_type=task_type, gt_h_file=gt_h_file,
                     exp_seed = exp_seed, covtype = "Matern5_2")

# save output
saveRDS(out, file=paste0("experiments", "/", exp_id, "_", "out.RData"))

# shut down EMEWS DB etc
if (!is.null(task_queue)) task_queue$close()
if (!is.null(pool)) pool$cancel()
if (db_started) eqsql$db_tools$stop_db(params$db_path)

