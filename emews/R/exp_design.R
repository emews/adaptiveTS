exp_design <- list(
  sim_budget = 3000, # Total simulation budget
  init_npar = 30,    # Number of initial design points
  nrep = 30,         # Number of replications
  grid_npar = 100,    # Number of candidate points for grid
  nTS_samp = 200,    # Number of Thompson samples
  p = 3,           # Dimension of the input space
  objective = "both", # one of (hosp, deaths, both)
  err_sig = 0.5,
  prop_sig = 0.3,
  ref = -2
)
