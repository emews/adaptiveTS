exp_design <- list(
  sim_budget = 5000, # Total simulation budget
  init_npar = 100,    # Number of initial design points
  nrep = 10,         # Number of replications
  grid_npar = 1000,    # Number of candidate points for grid
  nTS_samp = 100,    # Number of Thompson samples
  objective = "both", # one of (hosp, deaths, both)
  err_sig = 0.5,
  prop_sig = 0.3,
  ref = -2,
  param_names = c('infected.count', 
                  'susceptible.to.exposed.probability', 
                  'seasonality.multiplier',
                  'shielding.scaling', 
                  'isolate.infectivity.household',
                  'isolate.infectivity.nursinghome', 
                  'initial.exposure.tick',
                  'stay.at.home.probability', 
                  'stoe.behavioral.adjustment.probability')
)
