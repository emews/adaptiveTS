## surrogate

#' Title
#'
#' @param X 
#' @param Y 
#' @param GP_type 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
fitGP <- function(X, Y, GP_type, ...){
  if(GP_type == "homGP"){
    f <- mleHomGP(X[, 1:(ncol(X) - 1)], Y, ...)
  }
  if(GP_type == "hetGP"){
    f <- mleHetGP(X[, 1:(ncol(X) - 1)], Y, ...)
  }
  if(GP_type == "CRNGP"){
    f <- mleCRNGP(X, Y, ...)
  }
  return(f)
}

#' Title
#'
#' @param Ysim 
#' @param Ytrue 
#' @param err_sig 
#'
#' @return
#' @export
#'
#' @examples
loglik <- function(Ysim, Ytrue, err_sig){
  dnorm(Ysim, Ytrue, err_sig, log = T)
}



#' Title
#'
#' @param XSgrid 
#' @param model 
#' @param Ytrue 
#' @param err_sig 
#' @param check 
#'
#' @return
#' @export
#'
#' @examples
loglik_design <- function(model, newX, Ytrue, err_sig){
  
  # predict from fitted GP
  pred <- predict(model, newX, xprime = newX)
  
  ll <- loglik(pred$mean, Ytrue, err_sig)
  return(ll)
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
create_grid_CRNGP <- function(nparam, nrep, model, Ytrue, err_sig, prop_sig = 0.3){
  
  p <- ncol(model$X0)  # input dimension
  
  # lhs design to initialize
  s <- 1:nrep
  Xgrid_01 <- randomLHS(n = nparam, k = p)
  Xsgrid_01 <- cbind(Xgrid_01[rep(1:nparam, each = nrep), ], rep(s, nparam))
  
  # calculate importance weights
  w <- loglik_design(newX = Xsgrid_01, model = model, Ytrue = Ytrue, err_sig = err_sig)
  w_prime <-  w - max(w)
  w_norm <- exp(w_prime) / sum(exp(w_prime))
  
  # sample
  grid_ids <- sample(1:nrow(Xsgrid_01), nrow(Xsgrid_01), prob = w_norm, replace = T)
  grid_ids <- unique(grid_ids)
  
  # desnify
  xsgrid_new <- c()
  xsgrid_all <- c()
  xsgrid_id <- c()
  
  acc_prop <- c()
  
  k <- 1
  for (ii in grid_ids){
    
    if(nrow(xsgrid_new) < (nparam*nrep) || is.null(nrow(xsgrid_new))) {
      # old candidate
      xs_old <- Xsgrid_01[ii, ]
      # proposal
      xs_can <- rnorm(p, xs_old[1:p], prop_sig)
      # nrep_can <- sample.int(nrep, 1)
      xs_can_mat <- cbind(matrix(rep(xs_can, each = nrep), nrow = nrep), 1:nrep)
      
      # acceptance probability
      if((all(xs_can[1:p] < 1)) & ((all(xs_can[1:p] > 0)))){
        acc_prob <- loglik_design(xs_can_mat, model = f,
                                  Ytrue, .5) - w[ii]
        u <- log(runif(nrep))
        
        # print(length(u))
        # print(length(acc_prob))
        
        selected_ids <- which(u < acc_prob)
        if(length(selected_ids) > 0){
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
  
  return(list(xsgrid_new, xsgrid_all, xsgrid_id))
}


#' Title
#'
#' @param nparam 
#' @param nrep 
#' @param p 
#' @param model 
#' @param Ytrue 
#' @param err_sig 
#' @param prop_sig 
#'
#' @return
#' @export
#'
#' @examples
create_grid_GP <- function(nparam, nrep, p, model, Ytrue, err_sig, prop_sig = 0.3){
  
  p <- ncol(model$X0)  # input dimension
  
  # lhs design to initialize
  s <- 1:nrep
  Xgrid_01 <- randomLHS(n = nparam, k = p)
  Xsgrid_01 <- Xgrid_01[rep(1:nparam, each = nrep), ]
  
  # calculate importance weights
  w <- loglik_design(xsgrid = Xgrid_01, model = model, Ytrue = Ytrue, err_sig = err_sig)
  w <- rep(w, each = nrep)
  w_prime <-  w - max(w)
  w_norm <- exp(w_prime) / sum(exp(w_prime))
  
  # sample
  grid_ids <- sample(1:nrow(Xsgrid_01), nrow(Xsgrid_01), prob = w_norm, replace = T)
  grid_ids <- unique(grid_ids)
  
  # desnify
  xsgrid_new <- c()
  xsgrid_all <- c()
  xsgrid_id <- c()
  
  acc_prop <- c()
  
  k <- 1
  for (ii in grid_ids){
    
    if(nrow(xsgrid_new) < (nparam*nrep) || is.null(nrow(xsgrid_new))) {
      # old candidate
      xs_old <- Xsgrid_01[ii, ]
      # proposal
      xs_can <- rnorm(p, xs_old, prop_sig)
      # nrep_can <- sample.int(nrep, 1)
      xs_can_mat <- matrix(rep(xs_can, each = nrep), nrow = nrep)
      
      # acceptance probability
      if((all(xs_can[1:p] < 1)) & ((all(xs_can[1:p] > 0)))){
        loglik_new_design <- loglik_design(xs_can, model = f,
                                           Ytrue, err_sig)
        loglik_new_design <- rep(loglik_new_design, each = nrep)
        
        acc_prob <- loglik_new_design - w[ii]
        
        u <- log(runif(nrep))
        selected_ids <- which(u < acc_prob)
        
        if(length(selected_ids) > 0){
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
  
  return(list(xsgrid_new, xsgrid_all, xsgrid_id))
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
#'
#' @return
#' @export
#'
#' @examples
TSBatchBO <- function(init_npar,
                      nrep,
                      p, 
                      sim_budget,
                      grid_npar,
                      nTS_samp,
                      nTS_iter,
                      adaptive = TRUE,
                      ytrue,
                      GP_type,
                      sim_func,
                      exp_seed = NULL,
                      ...){
  
  
  
  ## =====================
  if(!is.null(exp_seed)) set.seed(exp_seed)
  
  ## =====================
  ## Fixed simulation settings
  nsteps <- 100
  beta_range <- c(0.2, 0.5)
  gamma_range <- c(0.1, 0.4)
  
  ## =====================
  ## initial design and simulations
  X_01 <- randomLHS(n = init_npar, k = p)
  s <- 1:nrep
  Xs_01 <- cbind(X_01[rep(1:init_npar, each = nrep), ], rep(s, init_npar))
  
  y <- rep(NA, init_npar*nrep)
  k <- 1
  for (ii in 1:init_npar){
    for (seed in s){
      y[k] <- sim_func(list(seed = seed, 
                            par = Xs_01[k, 1:p], 
                            param_bounds = rbind(beta_range, gamma_range), 
                            nsteps = nsteps,
                            ytrue_vec = ytrue))
      k <- k + 1
    }
  }
  
  ## =====================
  ## standardize output
  y_std <- scale(log(y))
  ycenter <- attr(y_std, "scaled:center")
  ysd <- attr(y_std, "scaled:scale")
  
  
  ## =====================
  ## initiate TS
  
  f <- fitGP(Xs_01, y_std, GP_type = GP_type)
  
  X_list <- list()
  y_list <- list()
  evaluated_ids <- 1:nrow(Xs_01)
  
  # ================================
  # fixed grid for non-adaptive case
  if(!adaptive) {
    Xgrid_01 <- randomLHS(n = grid_npar, k = p)
    Xsgrid_01 <- cbind(Xgrid_01[rep(1:grid_npar, each = nrep), ], 
                       rep(s, grid_npar))
    
    
    for (tt in 1:nTS_iter){
      out <- next_eval_CRN(model = f, 
                           Xsgrid = Xsgrid_01, 
                           evaluated_ids = evaluated_ids,
                           nTS_samp = nTS_samp,
                           adaptive = F)
      
      evaluated_ids <- c(evaluated_ids, out[[1]])
      xnew <- out[[2]]
      
      ## evaluate new simulations 
      ynew <- rep(NA, length(best_ids))
      
      if(!is.matrix(xnew)) xnew <- matrix(xnew, nrow = 1)
      
      for (ii in 1:length(best_ids)){
        ynew[ii] <- sim_func(list(seed = seed, 
                                  par = xnew[ii, 1:p], 
                                  param_bounds = rbind(beta_range, gamma_range), 
                                  nsteps = nsteps,
                                  ytrue_vec = ytrue))
      }
      X_list[[tt]] <- xnew
      y_list[[tt]] <- (log(ynew) - ycenter) / ysd  
      
      ## update surrogate
      f <- fitGP(rbind(Xs_01, X_list[[tt]]), c(y_std, Y_list[[tt]]), GP_type = GP_type)
    }
    
    return(list(X_list, Y_list))
  }
  
  
  
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
next_eval_CRN <- function(model, 
                          Xsgrid = NULL, 
                          evaluated_ids = NULL,
                          nTS_samp,
                          adaptive = TRUE,
                          grid_npar = NULL,
                          nrep = NULL,
                          Ytrue = NULL,
                          err_sig = NULL,
                          prop_sig = NULL,
                          ...){
  
  if(adaptive){
    out <- adaptive_CRN_TS(model = model,
                           grid_npar = grid_npar,
                           nrep = nrep,
                           Ytrue = Ytrue,
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


#' Title
#'
#' @param model 
#' @param Xsgrid 
#' @param evaluated_ids 
#' @param nTS_samp 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
fixed_CRN_TS <- function(model,
                         Xsgrid,
                         evaluated_ids,
                         nTS_samp, ...){
  
  
  Xfull <- rbind(cbind(model$X0, model$S0), Xsgrid)
  tTS <- simul(object = model, Xfull, ids = evaluated_ids, 
               nsim = nTS_samp, ...)
  
  best_ids <- apply(tTS, 2, function(x){
    m_tmp <- cbind(x, 1:length(x))
    m_tmp <- m_tmp[-evaluated_ids, ]
    ids <- which.min(m_tmp[, 1])
    native_ids <- m_tmp[ids, 2]
    return(native_ids)
  })
  
  best_ids <- unique(best_ids)
  
  return(list(best_ids, Xfull[best_ids, ]))
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
adaptive_CRN_TS <- function(model,
                            grid_npar,
                            nrep,
                            Ytrue,
                            err_sig,
                            prop_sig,
                            nTS_samp){
  
  ## create grid
  grid <- create_grid_CRNGP(grid_npar, 
                            nrep, 
                            model, 
                            Ytrue, 
                            err_sig, 
                            prop_sig)
  Xsgrid <- grid[[1]]
  Xsgrid_full <- grid[[2]]
  Xsgrid_id <- grid[[3]]
  
  ## predict
  pred <- predict(f, Xsgrid, xprime = Xsgrid)
  tTS <- MASS::mvrnorm(n = nTS_samp, 
                       mu = pred$mean, Sigma = 1/2 * (pred$cov + t(pred$cov)))
  
  best_ids <- apply(tTS, 1, which.min)
  best_ids <- unique(best_ids)
  
  return(Xsgrid[best_ids, ])
}
















