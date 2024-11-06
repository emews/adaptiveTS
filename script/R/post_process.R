## =================================================
## collect and format results

x_fixed_CRNGP_mat <- do.call(rbind, out_fixed_CRNGP[[1]])
x_adaptive_CRNGP_mat <- do.call(rbind, out_adaptive_CRNGP[[1]])
x_fixed_hetGP_mat <- do.call(rbind, out_fixed_hetGP[[1]])
x_adaptive_hetGP_mat <- do.call(rbind, out_adaptive_hetGP[[1]])

y_fixed_CRNGP_vec <- unlist(out_fixed_CRNGP[[2]])
y_adaptive_CRNGP_vec <- unlist(out_adaptive_CRNGP[[2]])
y_fixed_hetGP_vec <- unlist(out_fixed_hetGP[[2]])
y_adaptive_hetGP_vec <- unlist(out_adaptive_hetGP[[2]])

y_fixed_CRNGP_native <- unlist(out_fixed_CRNGP[[3]])
y_adaptive_CRNGP_native <- unlist(out_adaptive_CRNGP[[3]])

## run simulations

sim_length <- 101

y_fixed_CRNGP_traj <- matrix(NA, nrow = sim_length, ncol = nrow(x_fixed_CRNGP_mat))
for (ii in 1:nrow(x_fixed_CRNGP_mat)){
  o <- run_sim(x_fixed_mat[ii, ])
  y_fixed_traj[, ii] <- o$n_SI
}

y_adaptive_CRNGP_traj <- matrix(NA, nrow = sim_length, ncol = nrow(x_adaptive_CRNGP_mat))
for (ii in 1:nrow(x_adaptive_CRNGP_mat)){
  o <- run_sim(x_adaptive_CRNGP_mat[ii, ])
  y_adaptive_CRNGP_traj[, ii] <- o$n_SI
}

y_fixed_hetGP_traj <- matrix(NA, nrow = sim_length, ncol = nrow(x_fixed_hetGP_mat))
for (ii in 1:nrow(x_fixed_hetGP_mat)){
  o <- run_sim(x_fixed_hetGP_mat[ii, ])
  y_fixed_traj[, ii] <- o$n_SI
}

y_adaptive_hetGP_traj <- matrix(NA, nrow = sim_length, ncol = nrow(x_adaptive_hetGP_mat))
for (ii in 1:nrow(x_adaptive_hetGP_mat)){
  o <- run_sim(x_adaptive_hetGP_mat[ii, ])
  y_adaptive_hetGP_traj[, ii] <- o$n_SI
}

# run simulations and calculate sum squared errors

y_fixed_CRNGP_native <- y_adaptive_CRNGP_native <- y_fixed_hetGP_native <- y_adaptive_hetGP_native <- rep(NA, sim_budget)
for(ii in 1:sim_budget){
  y_fixed_CRNGP_native[ii] <-run_sim_err(x_fixed_CRNGP_mat[ii, ])
  y_adaptive_CRNGP_native[ii] <- run_sim_err(x_adaptive_CRNGP_mat[ii, ])
  y_fixed_hetGP_native[ii] <-run_sim_err(x_fixed_hetGP_mat[ii, ])
  y_adaptive_hetGP_native[ii] <- run_sim_err(x_adaptive_hetGP_mat[ii, ])
}


## ================================================
## plot

p_pairs <- pairs_plot(x_fixed_CRNGP_mat, 
                      x_adaptive_CRNGP_mat, 
                      x_fixed_hetGP_mat, 
                      x_adaptive_hetGP_mat)
p_pairs

for(ii in c(10, 30, 50, 80, 100)){
  print(plot_nbest_traj(nbest = ii, x_fixed_CRNGP_mat, x_adaptive_CRNGP_mat, 
                        y_fixed_CRNGP_vec, y_adaptive_CRNGP_vec, yobs = ytrue_vec))
  # print(err_hist(nbest = ii, y_fixed_native, y_adaptive_native))
}

for(ii in c(10, 30, 50, 80, 100)){
  print(plot_nbest_traj(nbest = ii, x_fixed_hetGP_mat, x_adaptive_hetGP_mat, 
                        y_fixed_hetGP_vec, y_adaptive_hetGP_vec, yobs = ytrue_vec))
  # print(err_hist(nbest = ii, y_fixed_native, y_adaptive_native))
}



nbests <- seq(20, sim_budget, length.out = 5)
par(mfrow = c(3, 2))
for(nbest in nbests){
  plot_pct_traj(nbest, y_fixed_CRNGP_native, y_adaptive_CRNGP_native)
}

par(mfrow = c(5, 2))
for(nbest in nbests){
  plot_pct_traj(nbest, y_fixed_hetGP_native, y_adaptive_hetGP_native)
}

par(mfrow = c(5, 2))
for(nbest in nbests){
  plot_pct_traj(nbest, y_adaptive_CRNGP_native, y_adaptive_hetGP_native,
                legend_labels = c("CRNGP-adaptive", "hetGP-adative"))
}


contour_compare(nbest = 50, 1, x_fixed_mat, x_adaptive_mat, y_fixed_vec, y_adaptive_vec, 
                d_dense[d_dense$seed == 1, ])
