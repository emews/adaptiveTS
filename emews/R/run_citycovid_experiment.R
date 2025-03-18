library(hetGP, lib.loc = "/lcrc/project/EMEWS/bebop/R-packages/")
library(scales)
library(data.table)
library(lhs)
library(ggplot2)
library(mvtnorm)
library(argparse)
library(reticulate)
library(EQ.SQL)

## =================================
## BO experiment settings 

# init_npar = 30,    # Number of initial design points
# nrep = 30,   

# hard-coding these for now, can make input if necessary
exp_design <- list(
  sim_budget = 3000, # Total simulation budget
  init_npar = 30,    # Number of initial design points
  nrep = 30,         # Number of replications
  grid_npar = 100,    # Number of candidate points for grid
  nTS_samp = 200,    # Number of Thompson samples
  p = 3,           # Dimension of the input space
  err_sig = 0.5,
  prop_sig = 0.3,
  ref = -2
)

## =================================
## EMEWS setup

run <- function(exp_id, params) {
  db_started <- FALSE
  pool <- NULL
  task_queue <- NULL

  tryCatch({
      eqsql <- init_eqsql(python_path = params$python_path)

      # eqsql$db_tools$stop_db(params$db_path)

      db_started <- FALSE
      pool <- NULL
      task_queue <- NULL

      eqsql <- init_eqsql(python_path = params$python_path)

      # eqsql$db_tools$start_db(params$db_path)
      # db_started <- TRUE

      task_queue <- init_task_queue(eqsql, params$db_host, params$db_user, params$db_port,
                                    params$db_name)

      if (!task_queue$are_queues_empty()) {
        print("WARNING: task input / output queues are not empty. Aborting run")
        task_queue$clear_queues()
      }
      pool_params <- eqsql$worker_pool$cfg_file_to_dict(params$pool_cfg_file)
      params$worker_pool_id <- 'adaptive_ts'
      exp_id <- 1
      if (params$pool_type == "local") {
          pool <- eqsql$worker_pool$start_local_pool(params$worker_pool_id, params$pool_launch_script,
                                                     exp_id, pool_params)
      } else {
          pool <- eqsql$worker_pool$start_scheduled_pool(params$worker_pool_id, params$pool_launch_script,
                                                         exp_id, pool_params, "pbs", gcx = NULL)
      }
      task_type <- params$task_type

      ## =================================
      ## run adaptive TS

      exp_seed <- 1 # Is this correct to set here?
      gt_h_file <- params$gt_h_file
      out <- runAdaptiveTS(exp_design, task_queue=task_queue, exp_id=exp_id, task_type=task_type, gt_h_file=gt_h_file,
                          exp_seed = exp_seed, covtype = "Matern5_2")

      # save output
      saveRDS(out, file=paste0(params$results_directory, "/", exp_id, "_", "out.RData"))
  }, error = function(e) {
      print(e)
      print(paste0("python error: ", reticulate::py_last_error()))
  }, finally = {
    # shut down EMEWS DB etc
    if (params$pool_type == "local") {
        print("Closing queue")
        if (!is.null(task_queue)) task_queue$close()
        print("Canceling pool")
        if (!is.null(pool)) pool$cancel()

    } else if (!is.null(task_queue)) {
        # pool$cancel seems to hang on improv with pbs
        print("Canceling pool")
        task_queue$stop_worker_pool(params$task_type)
        print("Closing queue")
        task_queue$close()
    }
        # if (db_started) eqsql$db_tools$stop_db(params$db_path)
  })
}

args <- commandArgs(trailingOnly = TRUE)
exp_id <- args[1]
params_file <- args[2]
params <- parse_yaml_cfg(params_file)
if (params$db_port == -1) {
    params$db_port = NULL
}

source(params$ts_r_file)

if (!file.exists(params$results_directory)) {
    dir.create(params$results_directory)
}

run(exp_id, params)
