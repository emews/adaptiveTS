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
  param_names = c('infected_count', 
                  'susceptible_to_exposed_probability', 
                  'seasonality_multiplier',
                  'shielding_scaling', 
                  'isolate_infectivity_household',
                  'isolate_infectivity_nursinghome', 
                  'initial_exposure_tick',
                  'stay_at_home_probability', 
                  'stoe_behavioral_adjustment_probability')
)
