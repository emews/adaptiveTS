library(odin)

## =====================================================================
## =====================================================================
## A basic SIR model in odin

stoch_SIR_odin <- odin::odin({
  ## this block of code looks like R, but they are actually
  ## parsed and compiled in C. 
  
  ## Core equations for transitions between compartments:
  
  update(S) <- S - n_SI
  update(I) <- I + n_SI - n_IR
  update(R) <- R + n_IR
  update(T) <- T + 1
  
  ## Individual probabilities of transition:
  p_SI <- 1 - exp(-beta * I / N) # S to I
  p_IR <- 1 - exp(-gamma) # I to R
  
  
  ## Draws from binomial distributions for numbers changing between
  ## compartments:
  n_SI <- if (n_SI_ini != -1 && T == 0) n_SI_ini else rbinom(S, p_SI)
  n_IR <- if (n_IR_ini != -1 && T == 0) n_IR_ini else rbinom(I, p_IR) 
  
  ## Total population size
  N <- S + I + R
  
  ## Initial states:
  initial(S) <- S_ini
  initial(I) <- I_ini
  initial(R) <- R_ini
  initial(T) <- 0
  
  ## User defined parameters - default in parentheses:
  S_ini <- user(1000)
  I_ini <- user(10)
  R_ini <- user(0)
  beta <- user(0.2)
  gamma <- user(0.1)
  n_SI_ini <- user(-1)
  n_IR_ini <- user(-1)
  
  output(n_SI) <- TRUE
  output(n_IR) <- TRUE
}, target = "c")



parse_sim_arg <- function(parlist){
  
  lapply(names(parlist), function(name) {
    assign(name, parlist[[name]], envir = sys.frame())
  })
  
}

## =====================================================================

#' Title
#'
#' @param seed 
#' @param params 
#' @param nsteps 
#' @param return_trajectory 
#' @param param_bounds 
#'
#' @return
#' @export
#'
#' @examples
run_sim <- function(par){
  
  # mandatory inputs
  # seed, par, param_bounds, nsteps = 100, return_trajectory = F
  
  # check if global random seed is present
  if (exists(".Random.seed", .GlobalEnv))
    oldseed <- .GlobalEnv$.Random.seed
  else
    oldseed <- NULL

  # parse input
  # parse_sim_arg(parlist)
  npar <- length(par) # must be 2 for now
  
  beta_range <- c(0.2, 0.5)
  gamma_range <- c(0.1, 0.4)
  nsteps <- 100
  sim_seed <- par[3]
  
  set.seed(sim_seed)
  beta <- beta_range[1] + par[1] * diff(range(beta_range))
  gamma <- gamma_range[1] + par[2] * diff(range(gamma_range))

  ## simulate
  model <- stoch_SIR_odin$new(S_ini = 10000, I_ini = 20,
                              beta = beta, gamma = gamma)
  out <- data.frame(model$run(step = 0:nsteps))

  # restore random seed
  if (!is.null(oldseed))
    .GlobalEnv$.Random.seed <- oldseed
  else
    rm(".Random.seed", envir = .GlobalEnv)

  return(out)
}

## =====================================================================

#' Title
#'
#' @param seed 
#' @param params 
#' @param nsteps 
#' @param ytrue_vec 
#' @param param_bounds 
#'
#' @return
#' @export
#'
#' @examples
run_sim_err <- function(par){
  
  # mandatory inputs
  # seed, params, param_bounds, nsteps, ytrue_vec
  
  # check if global random seed is present
  if (exists(".Random.seed", .GlobalEnv))
    oldseed <- .GlobalEnv$.Random.seed
  else
    oldseed <- NULL
  
  ## parse input
  # parse_sim_arg(parlist)
  npar <- length(par) # must be 2 for now
  
  beta_range <- c(0.2, 0.5)
  gamma_range <- c(0.1, 0.4)
  nsteps <- 100
  sim_seed <- par[3]
  
  set.seed(sim_seed)
  beta <- beta_range[1] + par[1] * diff(range(beta_range))
  gamma <- gamma_range[1] + par[2] * diff(range(gamma_range))
  
  ## simulate
  model <- stoch_SIR_odin$new(S_ini = 10000, I_ini = 20,
                              beta = beta, gamma = gamma)
  out <- data.frame(model$run(step = 0:nsteps))
  
  ## standardize the output
  
  # restore random seed
  if (!is.null(oldseed)) 
    .GlobalEnv$.Random.seed <- oldseed
  else
    rm(".Random.seed", envir = .GlobalEnv)
  
  return(list(out=sum((out$n_SI - ytrue_vec)^2), res=out))
}

