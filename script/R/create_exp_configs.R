
# BASE_DIR <- "/lcrc/project/EMEWS/afadikar/git/adaptiveTS/experiments/exp_gt_beta0.3_gamma0.8_seed50/"
BASE_DIR <- "/lcrc/project/EMEWS/afadikar/git/adaptiveTS/experiments/exp_gt_beta0.3_gamma0.8_seed50/"

if (!dir.exists(BASE_DIR)) {
  dir.create(BASE_DIR)
}

exp_seed <- 1:10
init_npar <- c(5, 10)
nrep <- c(10, 20)	
sim_budget <- c(700)	
grid_npar <- c(100, 200, 300)	
nTS_samp <- c(10, 20, 30)

d <- expand.grid(exp_seed, init_npar, nrep, sim_budget, grid_npar, nTS_samp)
colnames(d) <- c("exp_seed", "init_npar", "nrep", "sim_budget", "grid_npar", "nTS_samp")

d_eff <- d[d$init_npar * d$nrep < d$sim_budget, ]
d_eff$exp_id <- paste0("exp_", 1:nrow(d_eff))

d_eff <- d_eff[, c(7, 1:6)]
d_eff$exp_path <- paste0(BASE_DIR, d_eff$exp_id, "/")

UPF_DIR <- "/lcrc/project/EMEWS/afadikar/git/adaptiveTS/emews/data/upfs/"

write.csv(d_eff, file = paste0(UPF_DIR, "config_gt_beta0.3_gamma0.8_seed50.csv"), row.names = F)
