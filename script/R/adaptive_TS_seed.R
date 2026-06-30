

# Function to find matching row indices
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


# a function to create new grid at each iterations
#' Title
#'
#' @param nparam 
#' @param nrep 
#' @param model 
#' @param Ytrue 
#' @param err_sig 
#' @param prop_sig 
#'
#' @return
#' @export
#'
#' @examples
create_grid_CRNGP_seed <- function(evaluated_Xs, Xgrid_01, nrep){
  
  max_seed <- max(nrep, 
                  max(evaluated_Xs[, ncol(evaluated_Xs)]))
  
  # lhs design to initialize
  s <- 1:(max_seed+1)
  Xsgrid_01 <- cbind(Xgrid_01[rep(1:nrow(Xgrid_01), each = length(s)), ], 
                     rep(s, nrow(Xgrid_01)))
  
  return(Xsgrid_01)
}




parse_sim_arg <- function(parlist){
  
  
}



#' Title
#'
#' @param init_npar 
#' @param nrep 
#' @param p 
#' @param sim_budget 
#' @param grid_npar 
#' @param nTS_samp 
#' @param nTS_iter 
#' @param adaptive 
#' @param ytrue 
#' @param GP_type 
#' @param sim_func 
#' @param exp_seed 
#' @param ... 
#' @param ref 
#' @param err_sig 
#' @param prop_sig 
#'
#' @return
#' @export
#'
#' @examples
TSBatchBO_CRNGP_seed <- function(init_npar,
                                 nrep,
                                 p, 
                                 sim_budget,
                                 grid_npar,
                                 nTS_samp,
                                 nTS_iter = NULL,
                                 ytrue,
                                 GP_type = "CRNGP",
                                 sim_func,
                                 exp_seed = NULL,
                                 ...){
  
  
  
  ## =====================
  if(!is.null(exp_seed)) set.seed(exp_seed)
  
  ## =====================
  ## initial design and simulations
  X_01 <- lhs::randomLHS(n = init_npar, k = p)
  s <- 1:nrep
  Xs_01 <- cbind(X_01[rep(1:init_npar, each = nrep), ], rep(s, init_npar))
  
  simouts <- list() # store full simulations so we don't have to re-run later in the analysis
  y <- rep(NA, init_npar*nrep)
  for (ii in 1:(init_npar*nrep)){
    simout <- sim_func(Xs_01[ii,])
    y[ii] <- simout$out
    simouts[[ii]] <- simout$res
  }
  
  ## =====================
  ## standardize output
  y_std <- scale(log(y+1))
  ycenter <- attr(y_std, "scaled:center")
  ysd <- attr(y_std, "scaled:scale")
  
  
  ## =====================
  ## BO set up 
  Xs <- Xs_01
  Y <- y_std
  f <- fitGP(Xs, Y, GP_type = GP_type, known = list(beta0 = 0), ...)
  evaluated_Xs <- Xs
  
  no_of_sims <- length(y)
  X_list <- list()
  y_list <- list()
  simout_list <- list()
  ynative_list <- list()
  
  X_list[[1]] <- Xs_01
  y_list[[1]] <- y_std
  simout_list[[1]] <- simouts
  ynative_list[[1]] <- y
  
  # ================================
  # fixed grid for non-adaptive case

  # ================================
  # adaptive seed grid 
  if(TRUE){
    
    Xgrid_01 <- randomLHS(n = grid_npar, k = p)
    Xgrid_full_01 <- rbind(X_01, Xgrid_01)
    Xsgrid_01 <- cbind(Xgrid_01[rep(1:grid_npar, each = nrep), ], 
                       rep(s, grid_npar))
    Xsfull <- rbind(Xs_01, Xsgrid_01)
    
    # TS starts here
    tt <- 2
    while(no_of_sims < sim_budget){
      # for (tt in 1:nTS_iter){
      out <- next_eval_CRN_seed(model = f,
                                Xgrid_01 = Xgrid_full_01,
                                nTS_samp = nTS_samp,
                                evaluated_Xs = evaluated_Xs,
                                nrep = nrep)
      
      xnew <- out
      if(!is.matrix(xnew)) xnew <- matrix(xnew, nrow = 1)
      
      ## evaluate new simulations 
      ynew <- rep(NA, nrow(xnew))
      simouts <- list()
      
      for (ii in 1:nrow(xnew)){
        simout <- sim_func(xnew[ii, ])
        ynew[ii] <- simout$out
        simouts[[ii]] <- simout$res
      }
      X_list[[tt]] <- xnew
      y_list[[tt]] <- (log(ynew+1) - ycenter) / ysd
      ynative_list[[tt]] <- ynew
      simout_list[[tt]] <- simouts
      
      ## update surrogate
      Xs <- rbind(Xs, xnew)
      Y <- c(Y, y_list[[tt]])
      f <- fitGP(Xs, Y, GP_type = GP_type, 
                 known = list(beta0 = 0), ...)
      evaluated_Xs <- Xs
      
      no_of_sims <- no_of_sims + length(ynew)
      
      cat("iter = ", tt, "\n")
      
      tt <- tt + 1
    }
  }
  
  return(list("X_list"=X_list, "y_list"=y_list, "ynative_list"=ynative_list, "simout_list"=simout_list))
  
  # return(f)
  
}


#' Title
#'
#' @param model 
#' @param Xsgrid 
#' @param evaluated_ids 
#' @param nTS_samp 
#' @param adaptive 
#' @param grid_npar 
#' @param nrep 
#' @param Ytrue 
#' @param err_sig 
#' @param prop_sig 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
next_eval_CRN_seed <- function(model,
                               Xgrid_01,
                               nTS_samp,
                               evaluated_Xs,
                               nrep){
  
  if(TRUE){
    out <- adaptive_seed_CRN_TS(model = model, 
                                evaluated_Xs = evaluated_Xs,
                                Xgrid_01 = Xgrid_01,
                                nTS_samp = nTS_samp,
                                nrep = nrep)
    
  } 
  return(out)
}


#' Title
#'
#' @param model 
#' @param grid_npar 
#' @param nrep 
#' @param Ytrue 
#' @param err_sig 
#' @param prop_sig 
#' @param nTS_samp 
#'
#' @return
#' @export
#'
#' @examples
adaptive_seed_CRN_TS <- function(model,
                                 evaluated_Xs,
                                 Xgrid_01,
                                 nrep,
                                 nTS_samp){
  
  ## create grid
  eff_grid <- create_grid_CRNGP_seed(evaluated_Xs = evaluated_Xs,
                                     Xgrid_01 = Xgrid_01,
                                     nrep = nrep)
  
  Xsgrid <- eff_grid
  
  ## remove already evaluted points
  evaluated_ids <- find_row_indices(Xsgrid, evaluated_Xs)
  # eff_grid <- Xsgrid[-evaluated_ids, ]
  
  ## predict
  # pred <- predict(model, eff_grid, xprime = eff_grid)
  # tTS <- MASS::mvrnorm(n = nTS_samp,
  #                      mu = pred$mean,
  #                      Sigma = 1/2 * (pred$cov + t(pred$cov)))
  tTS <- simul(object = model, Xsgrid, ids = evaluated_ids,
               nsim = nTS_samp, check = F)
  
  # best_ids <- apply(tTS, 1, which.min)
  best_ids <- apply(tTS, 2, function(x){
    m_tmp <- cbind(x, 1:length(x))
    m_tmp <- m_tmp[-evaluated_ids, ]
    ids <- which.min(m_tmp[, 1])
    native_ids <- m_tmp[ids, 2]
    return(native_ids)
  })
  best_ids <- sort(unique(best_ids))
  
  return(eff_grid[best_ids, ])
}


