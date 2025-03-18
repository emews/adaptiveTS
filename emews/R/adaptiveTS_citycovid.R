# Combined Adaptive Thompson Sampling (TS) Framework
# ============================================================

# Load required libraries
library(hetGP)
library(scales)
library(data.table)
library(lhs)
library(ggplot2)
library(mvtnorm)
library(argparse)


# -----------------------------------------------------------
# Core GP Modeling Functions
# -----------------------------------------------------------

#' Fit a Gaussian Process model of specified type
#'
#' @param X Design matrix with input variables and optional seed column
#' @param Y Response values
#' @param GP_type Type of GP model: "homGP", "hetGP", or "CRNGP"
#' @param ... Additional parameters to pass to the specific GP fitting function
#'
#' @return Fitted GP model
#' @export
fitGP <- function(X, Y, GP_type, ...) {
  if (GP_type == "homGP") {
    f <- mleHomGP(X[, 1:(ncol(X) - 1)], Y, ...)
  } else if (GP_type == "hetGP") {
    f <- mleHetGP(X[, 1:(ncol(X) - 1)], Y, ...)
  } else if (GP_type == "CRNGP") {
    f <- mleCRNGP(X, Y, ...)
  } else {
    stop("Unsupported GP_type: ", GP_type)
  }
  return(f)
}

#' Calculate log-likelihood for comparing simulated outputs to true values
#'
#' @param ysim Simulated or predicted values
#' @param ytrue True values or reference values
#' @param err_sig Error standard deviation
#'
#' @return Log-likelihood values
#' @export
loglik <- function(ysim, ytrue, err_sig) {
  dnorm(ysim, ytrue, err_sig, log = TRUE)
}

# -----------------------------------------------------------
# Grid Creation Functions
# -----------------------------------------------------------

#' Create an adaptive grid for CRN GP models
#'
#' @param nparam Number of parameter combinations
#' @param nrep Number of replications for each parameter set
#' @param model Fitted GP model
#' @param ref Reference values for comparison
#' @param err_sig Error standard deviation
#' @param prop_sig Proposal standard deviation for the MCMC step
#'
#' @return List with grid matrices
#' @export
create_grid_CRNGP <- function(nparam, nrep, model, ref, err_sig, prop_sig = 0.3) {
  p <- ncol(model$X0)  # input dimension
  
  # LHS design to initialize
  s <- 1:nrep
  Xgrid_01 <- randomLHS(n = nparam, k = p)
  Xsgrid_01 <- cbind(Xgrid_01[rep(1:nparam, each = nrep), ], rep(s, nparam))
  
  # Calculate importance weights
  w <- loglik_design(newX = Xsgrid_01, model = model, ytrue = ref, err_sig = err_sig)
  w_prime <- w - max(w)
  w_norm <- exp(w_prime) / sum(exp(w_prime))
  
  # Sample
  grid_ids <- sample(1:nrow(Xsgrid_01), nrow(Xsgrid_01), prob = w_norm, replace = TRUE)
  grid_ids <- unique(grid_ids)
  
  # Densify
  xsgrid_new <- c()
  xsgrid_all <- c()
  xsgrid_id <- c()
  acc_prop <- c()
  
  k <- 1
  for (ii in grid_ids) {
    if (nrow(xsgrid_new) < (nparam * nrep) || is.null(nrow(xsgrid_new))) {
      # Old candidate
      xs_old <- Xsgrid_01[ii, ]
      # Proposal
      xs_can <- rnorm(p, xs_old[1:p], prop_sig)
      # Generate replications with seed IDs
      xs_can_mat <- cbind(matrix(rep(xs_can, each = nrep), nrow = nrep), 1:nrep)
      
      # Acceptance probability
      if ((all(xs_can[1:p] < 1)) & (all(xs_can[1:p] > 0))) {
        acc_prob <- loglik_design(xs_can_mat, model = model, ref, err_sig) - w[ii]
        u <- log(runif(nrep))
        
        selected_ids <- which(u < acc_prob)
        if (length(selected_ids) > 0) {
          selected_id <- sort(sample(selected_ids, ceiling(length(selected_ids)/2)))
          xsgrid_new <- rbind(xsgrid_new, xs_can_mat[selected_id, ])
          xsgrid_all <- rbind(xsgrid_all, xs_can_mat)
          xsgrid_id <- c(xsgrid_id, selected_id + (k-1)*20)
        }
        k <- k + 1
        acc_prop <- c(acc_prop, 1)
      } else {
        acc_prop <- c(acc_prop, 0)
      }
    }
  }

  # print(paste0("xsgrid_new: ", xsgrid_new))
  
  return(list(xsgrid_new, xsgrid_all, xsgrid_id))
}

# -----------------------------------------------------------
# Log-Likelihood Evaluation Functions
# -----------------------------------------------------------

#' Evaluate log-likelihood for design
#'
#' @param model Fitted GP model
#' @param newX New input points
#' @param ytrue Reference values
#' @param err_sig Error standard deviation
#'
#' @return Log-likelihood values
#' @export
loglik_design <- function(model, newX, ytrue, err_sig) {
  # Predict from fitted GP
  pred <- predict(model, newX, xprime = newX)
  
  ll <- loglik(pred$mean, ytrue, err_sig)
  return(ll)
}

# -----------------------------------------------------------
# Next Point Selection Functions
# -----------------------------------------------------------

#' Select next points to evaluate for CRN TS
#'
#' @param model Fitted GP model
#' @param Xsgrid Optional grid of candidate points
#' @param evaluated_ids Optional indices of already evaluated points
#' @param nTS_samp Number of Thompson samples
#' @param adaptive Whether to use adaptive grid
#' @param grid_npar Number of parameter combinations for adaptive grid
#' @param nrep Number of replications
#' @param ref Reference values
#' @param err_sig Error standard deviation
#' @param prop_sig Proposal standard deviation
#' @param ... Additional parameters
#'
#' @return New points to evaluate
#' @export
next_eval_CRN <- function(model, 
                          Xsgrid = NULL, 
                          evaluated_ids = NULL,
                          nTS_samp,
                          adaptive = TRUE,
                          grid_npar = NULL,
                          nrep = NULL,
                          ref = NULL,
                          err_sig = NULL,
                          prop_sig = NULL,
                          ...) {
  
  if (adaptive) {
    out <- adaptive_CRN_TS(model = model,
                           grid_npar = grid_npar,
                           nrep = nrep,
                           ref = ref,
                           err_sig = err_sig,
                           prop_sig = prop_sig,
                           nTS_samp = nTS_samp, ...)
  } else {
    out <- fixed_CRN_TS(model = model,
                        Xsgrid = Xsgrid,
                        evaluated_ids = evaluated_ids,
                        nTS_samp = nTS_samp, ...)
  }
  
  return(out)
}


# -----------------------------------------------------------
# Thompson Sampling Implementation Functions
# -----------------------------------------------------------

#' Adaptive Thompson sampling for CRN models
#'
#' @param model Fitted GP model
#' @param grid_npar Number of parameter combinations
#' @param nrep Number of replications
#' @param ref Reference values
#' @param err_sig Error standard deviation
#' @param prop_sig Proposal standard deviation
#' @param nTS_samp Number of Thompson samples
#'
#' @return New points to evaluate
#' @export
adaptive_CRN_TS <- function(model,
                            grid_npar,
                            nrep,
                            ref,
                            err_sig,
                            prop_sig,
                            nTS_samp) {
  
  # Create grid
  grid <- create_grid_CRNGP(grid_npar, 
                            nrep, 
                            model, 
                            ref, 
                            err_sig, 
                            prop_sig)
  Xsgrid <- grid[[1]]
  Xsgrid_full <- grid[[2]]
  Xsgrid_id <- grid[[3]]
  
  # Predict
  pred <- predict(model, Xsgrid, xprime = Xsgrid)
  tTS <- MASS::mvrnorm(n = nTS_samp, 
                       mu = pred$mean, Sigma = 1/2 * (pred$cov + t(pred$cov)))
  
  best_ids <- apply(tTS, 1, which.min)
  best_ids <- unique(best_ids)
  
  return(Xsgrid[best_ids, ])
}


# -----------------------------------------------------------
# Utility Functions
# -----------------------------------------------------------

#' Find indices of matching rows in a larger matrix
#'
#' @param bigger_matrix Larger matrix to search in
#' @param smaller_matrix Smaller matrix to find
#'
#' @return Vector of indices
#' @export
find_row_indices <- function(bigger_matrix, smaller_matrix) {
  match_indices <- c()
  
  # Loop through each row in the smaller matrix
  for (i in 1:nrow(smaller_matrix)) {
    row_match <- apply(bigger_matrix, 1, 
                       function(row) all(row == smaller_matrix[i, ]))
    match_indices <- c(match_indices, which(row_match))
  }
  
  return(match_indices)
}

scale_params <- function(x, bounds){
  return(bounds[1] + x*(bounds[2]-bounds[1]))
}

#' Format parameters to EMEWS payload
#'
#' @param params vector of parameters to update payload

#' @return json string in correct format
generate_payload <- function(params, 
                             bounds=list(c(.01, .15), 
                                         c(0.1, 1),
                                         c(0.01, 0.5))) {
  # fixed
  fixed_params <- list(
    "infected.count" = 443.6,
    "seasonality.multiplier" = 0.211,
    "isolate.infectivity.household" = 0.471,
    "shielding.scaling" = 0.967,
    "isolate.infectivity.nursinghome" = 0.549,
    "initial.exposure.tick" = 72 #day 3, right?
  )
  
  # scaled inputs (must be in this order)
  var_params <- list(
    "susceptible.to.exposed.probability" = scale_params(params[1], bounds[[1]]),
    "stay.at.home.probability" = scale_params(params[2], bounds[[2]]),
    "stoe.behavioral.adjustment.probability" = scale_params(params[3], bounds[[3]]),
    "seed" = params[4]
  )
  
  # combine
  all_params <- c(fixed_params, var_params)
  
  # json
  json_string <- jsonlite::toJSON(all_params, auto_unbox = TRUE)
  
  return(json_string)
}

# -----------------------------------------------------------
# EMEWS submission
# -----------------------------------------------------------


#' Evaluate design points using EMEWS
#'
#' @param design_points Matrix of design points to evaluate
#' @param task_queue EMEWS task queue
#' @param exp_id Experiment ID
#' @param task_type Task type
#' @param iter_name Iteration name/identifier
#' @param ytrue True value for file-based loss calculations
#'
#' @return List containing outputs, simulation details, etc.
#' @export
submit_emews <- function(design_points, gt_h_file, task_queue, exp_id, task_type) {
  # Submit all tasks to EMEWS
  fts <- apply(design_points, 1, function(a) {
    payload <- generate_payload(a)
    submission <- task_queue$submit_task(exp_id, task_type, payload)
    submission[[2]]
  })
  
  # Collect results as they complete
  results <- as_completed(task_queue, fts, function(ft) {
    list(
      eq_task_id = ft$eq_task_id,
      data = ft$result()[[2]]
    )
  })
  
  # Process results
  eq_ids <- sapply(results$f_results, function(x) x$eq_task_id)
  eq_ids <- unlist(eq_ids)
  outfiles <- sapply(results$f_results[order(eq_ids)], function(x) x$data)
  # print(outfiles)
  
  # Calculate loss/output values from result files
  y <- unlist(lapply(outfiles, obj_h, gt_h_file = gt_h_file))
  # print(y)
  
  return(y)
}



#' Objective on hospitalization
#' # output_file <- "/lcrc/project/EMEWS/improv/ncollier/repos/snl_brave_citycovid/experiments/stein_mcmc_01212025_1.0/instances/instance_0/output/counts_r0_0.csv"
# gt_h_file <- "/lcrc/project/EMEWS/afadikar/git/CityCOVID_DA/data/dt.chicago.hosp.csv"
# gt_d_file <- "/lcrc/project/EMEWS/afadikar/git/CityCOVID_DA/data/dt.chicago.deaths.csv"
#'
#' @param output_file path to the CityCOVID count file
#' @param gt_h_file path to the chicago ground truth hospitalization
#' @param start_date beginning date of calibration period
#' @param end_date end date of calibration period
#'
#' @return sum of squared difference between simulated and observed trajectory

obj_h <- function(output_file, gt_h_file, 
                  start_date = as.Date('2020-03-17'), 
                  end_date = as.Date("2020-06-13")){
  
  ## process simulation
  sim_df <- data.table::fread(output_file,
                              select = c("tick", "hosp_r_count", "hosp_icu_r_count",
                                         "hosp_d_count", "hosp_icu_d_count",
                                         "icu_r_count", "icu_d_count"))
  sim_df[, date := as.IDate("2020-03-16") + (tick / 24)]
  sim_df[, total_hosp_sim := rowSums(.SD), .SDcols = patterns("count$")]
  
  ## process ground truth
  gt_df <- data.table::fread(gt_h_file)
  gt_df <- gt_df[(date >= start_date) & (date <= end_date)]
  gt_df <- gt_df[!is.na(tot.hosp)]
  
  ## merge
  df_all <- merge(sim_df, gt_df, by = "date", all.y = TRUE)
  
  return(sum((df_all$tot.hosp - df_all$total_hosp_sim)^2))
}

# -----------------------------------------------------------
# Main Optimization Functions
# -----------------------------------------------------------

#' Adaptive CRNGP TS
#'
#' @param exp_design Named list of experimental settings
#' @param ref Reference values
#' @param err_sig Error standard deviation
#' @param prop_sig Proposal standard deviation
#' @param task_queue EMEWS task queue (for EMEWS evaluation)
#' @param exp_id Experiment ID (for EMEWS evaluation)
#' @param task_type Task type (for EMEWS evaluation)
#' @param gt_h_file Filepath for data we are calibrating to
#' @param exp_seed Seed for reproducibility
#' @param ... Additional parameters
#'
#' @return List with optimization results
#' @export
runAdaptiveTS <- function(exp_design,
                          task_queue = NULL,
                          exp_id = NULL,
                          task_type = NULL,
                          gt_h_file = NULL,
                          exp_seed = NULL,
                          ...) {
  
  # Set seed if provided
  if (!is.null(exp_seed)) set.seed(exp_seed)
  sim_budget <- exp_design$sim_budget
  init_npar <- exp_design$init_npar
  nrep <- exp_design$nrep
  nTS_samp <- exp_design$nTS_samp
  grid_npar <- exp_design$grid_npar
  p <- exp_design$p
  prop_sig <- exp_design$prop_sig
  err_sig <- exp_design$err_sig
  ref <- exp_design$ref
  
  # Create initial design
  X_01 <- randomLHS(n = init_npar, k = p)
  s <- 1:nrep
  Xs_01 <- cbind(X_01[rep(1:init_npar, each = nrep), ], rep(s, init_npar))
  
  # Evaluate initial design
  print("Submitting Initial Design")
  y <- submit_emews(Xs_01, gt_h_file, task_queue, exp_id, task_type)
  # print("Finished Initial Design Evaluation")
  
  # Standardize output
  y_std <- scale(log(y))
  ycenter <- attr(y_std, "scaled:center")
  ysd <- attr(y_std, "scaled:scale")

  # Initialize BO
  Xs <- Xs_01
  Y <- y_std
  f <- fitGP(Xs, Y, GP_type = "CRNGP", known = list(beta0 = 0), ...)
  
  # Track iterations
  no_of_sims <- length(y)
  X_list <- list()
  y_list <- list()
  ynative_list <- list()

  X_list[[1]] <- Xs_01
  y_list[[1]] <- y_std
  ynative_list[[1]] <- y

  
  # ================================
  # adaptive grid 
  
  # TS starts here
  tt <- 2
  cat("Starting Iterative Evaluations\n")
  while(no_of_sims < sim_budget){
    # for (tt in 1:nTS_iter){
    print(f)
    cat("Num Sims: ", no_of_sims, "\n")
    out <- next_eval_CRN(model = f,
                         grid_npar = grid_npar,
                         nrep = nrep,
                         ref = ref,
                         err_sig = err_sig,
                         prop_sig = prop_sig,
                         nTS_samp = nTS_samp,
                         adaptive = T)
    
    xnew <- out
    # cat("xnew: ", xnew, "\n")
    flush.console()
    if(!is.matrix(xnew)) xnew <- matrix(xnew, nrow = 1)
    # print(paste0("xnew_m: ", xnew))
    
    ## evaluate new simulations 
    ynew <- submit_emews(xnew, gt_h_file, task_queue, exp_id, task_type)
    X_list[[tt]] <- xnew
    y_list[[tt]] <- (log(ynew) - ycenter) / ysd
    ynative_list[[tt]] <- ynew

    ## update surrogate
    Xs <- rbind(Xs, xnew)
    Y <- c(Y, y_list[[tt]])
    f <- fitGP(Xs, Y, GP_type = "CRNGP", 
               known = list(beta0 = 0), ...)
    
    no_of_sims <- no_of_sims + length(ynew)
    
    cat("iter = ", tt, "\n")
    
    tt <- tt + 1
  }
  
  return(list("X_list"=X_list, 
              "y_list"=y_list, 
              "ynative_list"=ynative_list,
              "final_model"=f,
              "standardization" = list("center" = ycenter, "scale" = ysd)))

}
