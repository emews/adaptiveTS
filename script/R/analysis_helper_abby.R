library(patchwork)
library(ggplot2)

collect_results <- function(out_list){
  x_mat <- do.call(rbind, out_list$X_list)
  y_vec <- unlist(out_list$y_list)
  y_native <- unlist(out_list$ynative_list)
  
  simout_flat <- do.call(c,out_list$simout_list)
  simout_mat <- do.call(cbind, lapply(simout_flat, function(df) df$n_SI))
  
  return(list(x_mat=x_mat, 
              y_vec=y_vec, 
              y_native=y_native,
              y_rmse=sqrt(y_native/nrow(simout_mat)),
              simout_mat=simout_mat,
              iter_idx=unlist(mapply(rep, seq_along(out_list$X_list), sapply(out_list$X_list, nrow)))))
}

collect_seed_results <- function(out_list){
  names(out_list) <- c('X_list', 'y_list', 'ynative_list')
  x_mat <- do.call(rbind, out_list$X_list)
  y_vec <- unlist(out_list$y_list)
  y_native <- unlist(out_list$ynative_list)
  
  # your dumbass forgot to save the full simulation outputs so maybe do that now lol
  simouts <- lapply(1:nrow(x_mat), function(ii) run_sim(x_mat[ii, ]))
  simout_mat <- do.call(cbind, lapply(simouts, function(df) df$n_SI))
  
  return(list(x_mat=x_mat, 
              y_vec=y_vec, 
              y_native=y_native,
              y_rmse=sqrt(y_native/nrow(simout_mat)),
              simout_mat=simout_mat,
              iter_idx=unlist(mapply(rep, seq_along(out_list$X_list), sapply(out_list$X_list, nrow)))))
}

load5outputs <- function(experiment_id, exp_path, config){
  outfile <- sprintf("%s/exp_%s_out.rds", exp_path, experiment_id)
  # only want to load the relevant elements, not the entire workspace
  results_list <- readRDS(outfile)
  
  res_fixed_CRNGP <- collect_results(results_list$out_fixed_CRNGP)
  res_adaptive_CRNGP <- collect_results(results_list$out_adaptive_CRNGP)
  res_fixed_hetGP <- collect_results(results_list$out_fixed_hetGP)
  res_adaptive_hetGP <- collect_results(results_list$out_adaptive_hetGP)
  res_adaptive_CRNGP_seed <- collect_seed_results(results_list$out_adaptive_CRNGP_seed)
  
  return(list(fixed_crn=res_fixed_CRNGP, 
              adaptive_crn=res_adaptive_CRNGP, 
              fixed_het=res_fixed_hetGP, 
              adaptive_het=res_adaptive_hetGP,
              fixed_crn_seed=res_adaptive_CRNGP_seed,
              config=config[exp_id==experiment_id,]))
}

load_outputs <- function(experiment_id, exp_path, config){
  outfile <- sprintf("%s/exp_%s_out.rds", exp_path, experiment_id)
  # only want to load the relevant elements, not the entire workspace
  results_list <- readRDS(outfile)
  
  res_fixed_CRNGP <- collect_results(out_fixed_CRNGP)
  res_adaptive_CRNGP <- collect_results(out_adaptive_CRNGP)
  res_fixed_hetGP <- collect_results(out_fixed_hetGP)
  res_adaptive_hetGP <- collect_results(out_adaptive_hetGP)
  res_adaptive_CRNGP_seed <- collect_results(out_adaptive_CRNGP_seed)
  
  return(list(fixed_crn=res_fixed_CRNGP, 
              adaptive_crn=res_adaptive_CRNGP, 
              fixed_het=res_fixed_hetGP, 
              adaptive_het=res_adaptive_hetGP,
              fixed_crn_seed=res_adaptive_CRNGP_seed,
              config=config[exp_id==experiment_id,]))
}


threshold_counts <- function(out_list, sim_budget=200, nout=5, thresholds=c(15, 20, 25, 30)) {
  outs <- names(out_list)[1:nout]
  all_results <- list()
  k = 1
  for (i in seq_along(outs)) {
    out_name <- outs[i]
    out <- out_list[[out_name]]
    for (j in seq_along(thresholds)) {
      t <- thresholds[j]
      cnt <- sum(out$y_rmse[1:sim_budget]<t)
      res <- list(method=out_name, threshold=t, num_trajectories=cnt)
      res <- cbind(data.frame(res), out_list$config)
      res$sim_budget <- sim_budget # overwrite
      all_results[[k]] <- res
      k = k+1
    }
  }
  return(do.call(rbind, all_results))
}

# load_outputs <- function(experiment_id, exp_path, config){
#   outfile <- sprintf("%s/exp_%s_out.RData", exp_path, experiment_id)
#   # only want to load the relevant elements, not the entire workspace
#   temp_env <- new.env()
#   load(outfile, envir = temp_env)
#   out_fixed_CRNGP <- temp_env$out_fixed_CRNGP
#   out_adaptive_CRNGP <- temp_env$out_adaptive_CRNGP
#   out_fixed_hetGP <- temp_env$out_fixed_hetGP
#   out_adaptive_hetGP <- temp_env$out_adaptive_hetGP
#   out_adaptive_CRNGP_seed <- temp_env$out_adaptive_CRNGP_seed
#   
#   res_fixed_CRNGP <- collect_results(out_fixed_CRNGP)
#   res_adaptive_CRNGP <- collect_results(out_adaptive_CRNGP)
#   res_fixed_hetGP <- collect_results(out_fixed_hetGP)
#   res_adaptive_hetGP <- collect_results(out_adaptive_hetGP)
#   res_adaptive_CRNGP_seed <- collect_results(out_adaptive_CRNGP_seed)
#   
#   return(list(fixed_crn=res_fixed_CRNGP, 
#               adaptive_crn=res_adaptive_CRNGP, 
#               fixed_het=res_fixed_hetGP, 
#               adaptive_het=res_adaptive_hetGP,
#               fixed_crn_seed=res_adaptive_CRNGP_seed,
#               config=config[exp_id==experiment_id,]))
# }

# threshold_counts <- function(out_list, nout=4, thresholds=c(15, 20, 25, 30)) {
#   outs <- names(out_list)[1:nout]
#   all_results <- list()
#   k = 1
#   for (i in seq_along(outs)) {
#     out_name <- outs[i]
#     out <- out_list[[out_name]]
#     for (j in seq_along(thresholds)) {
#       t <- thresholds[j]
#       cnt <- sum(out$y_rmse<t)
#       res <- list(method=out_name, threshold=t, num_trajectories=cnt)
#       all_results[[k]] <- cbind(data.frame(res), out_list$config)
#       k = k+1
#     }
#   }
#   return(do.call(rbind, all_results))
# }

checkpointing_speed <- function(out_list,nout=4, n_trajectories=c(10, 25, 50), thresholds=c(15, 20, 25, 30)) {
  outs <- names(out_list)[1:nout]
  all_results <- list()
  k = 1
  for (i in seq_along(outs)) {
    out_name <- outs[i]
    out <- out_list[[out_name]]
    for (t in thresholds) {
      for (n in n_trajectories){
        cnts <- which(cumsum(out$y_rmse<t) >= n)
        if (length(cnts) == 0){
          cnt <- -1
        }
        else {
          cnt <- min(cnts)
        }
        res <- list(method=out_name, threshold=t, n_trajectories=n, checkpoint=cnt)
        all_results[[k]] <- cbind(data.frame(res), out_list$config)
        k = k+1
      }
    }
  }
  return(do.call(rbind, all_results))
}

threshold_budget_analysis <- function(out_list, nout=4, thresholds=c(15, 20, 25, 30), show_plots=FALSE) {
  outs <- names(out_list)[1:nout]
  all_results <- list()
  budget_results <- list()
  k = 1
  for (i in seq_along(outs)) {
    out_name <- outs[i]
    out <- out_list[[out_name]]
    for (t in thresholds) {
      within_threshold <- cumsum(out$y_rmse<t)
      budget_df <- data.frame(time=1:length(within_threshold), trajectories_found=within_threshold)
      budget_df$method <- out_name
      budget_df$exp_id <- out_list$config$exp_id
      budget_df$threshold=t
      auc <- auc_trapezoids(within_threshold)
      res <- list(method=out_name, threshold=t, auc=auc)
      all_results[[k]] <- cbind(data.frame(res), out_list$config)
      budget_results[[k]] <- budget_df
      k = k+1
    }
  }
  budget_df <- do.call(rbind, budget_results)
  if(show_plots){
    print(ggplot(budget_df) + geom_line(aes(x=time, y=trajectories_found, color=method, group=method)) + facet_wrap(~threshold, ncol=4) + theme_minimal())
  }
  return(list(agg_auc=do.call(rbind, all_results),
              budget_df=budget_df))
}

analyze_experiment <- function(exp_id, exp_path, config, nout=5, thresholds=c(15, 20, 25, 30), n_trajectories=c(10, 25, 50)) {
  out_list <- load_outputs(exp_id, exp_path, config)
  # threshold analysis
  df_threshold <- threshold_counts(out_list, thresholds=thresholds, nout=nout)
  df_checkpoint <- checkpointing_speed(out_list, nout=nout, n_trajectories=n_trajectories, thresholds=thresholds)
  budg <- threshold_budget_analysis(out_list, nout=nout, thresholds=thresholds)
  return(list(threshold_count = df_threshold,
              checkpoint_reached = df_checkpoint,
              threshold_budget = budg$agg_auc,
              budget_traj = budg$budget_df))
}



nbest_mean <- function(nbest, out){
  y <- out$y_rmse
  ybest <- y[order(y)][1:nbest]
  return(mean(ybest))
}


auc_trapezoids<- function(y){
  x <- seq_along(y)
  # want zero to be reference
  auc <- sum((y[-1] + y[-length(y)]) * diff(x) / 2)
  return(auc)
}

nbest_pct <- function(nbest, y, ymax, nerror_ticks = 20){
    ybest <- y[order(y)][1:nbest]
    error_ticks <- seq(0, ymax, length.out=nerror_ticks)
    pctraj <- rep(NA, nerror_ticks)
    for (i in 1:nerror_ticks){
      pctraj[i] <- mean(ybest < error_ticks[i])
    }
    auc <- auc_trapezoids(pctraj)
    return(list(pctraj=pctraj,
                auc=auc,
                auc_pct=auc/nerror_ticks)
           )
}

#nbest_so_far <- function(nbest, y)

plot_thresholded_traj <- function(t, out, yobs){
  df_true <- data.frame(time= (1:length(yobs)) - 1, yobs = yobs) 
  ix <- which(out$y_rmse < t)
  if(length(ix)==0){
    ggplot() +
      geom_line(data = df_true, mapping = aes(x = time, y = yobs), color = "black", linewidth = 1) +
      annotate("text", x = max(df_true$time) * 0.8, y = max(df_true$yobs) * 0.8, 
               label = paste("# traj:", length(ix))) +
      theme_minimal()
  }
  else{
    df_traj <- as.data.frame(out$simout_mat[, ix])
    df_traj$time <- (1:nrow(df_traj)) - 1
    df_traj <- pivot_longer(df_traj, 
                   cols = -time, 
                   names_to = "traj", 
                   values_to = "new_infections")
    
    ggplot() +
      geom_line(data=df_traj, aes(x=time, y=new_infections, group=traj), color="#F8766D", alpha=0.5) +
      geom_line(data = df_true, mapping = aes(x = time, y = yobs), color = "black", linewidth = 1) +
      annotate("text", x = max(df_traj$time) * 0.8, y = max(df_traj$new_infections) * 0.8, 
               label = paste("# traj:", length(ix))) +
      theme_minimal()
  }
}

plot_thresholds <- function(out_list, yobs, thresholds=c(15, 20, 25, 30)){
  outs <- names(out_list)[1:4]  # Assuming out_list is a named list
  config_info <- out_list$config
  main_title <- paste("Experiment", config_info$exp_id)
  subtitle <- paste(
    "Random seed =", config_info$exp_seed,
    "| Initial # parameters =", config_info$init_npar,
    "| Initial # reps =", config_info$nrep,
    "| Budget =", config_info$sim_budget,
    "| Grid size =", config_info$grid_npar,
    "| TS sample size =", config_info$nTS_samp
  )
  
  # Create a list to store all plots
  plots <- list()
  
  # Loop over each combination of thresholds and outs
  for (i in seq_along(outs)) {
    out_name <- outs[i]
    for (j in seq_along(thresholds)) {
      t <- thresholds[j]
      
      # Generate the plot
      p <- plot_thresholded_traj(t, out_list[[out_name]], ytrue) +
        ggtitle(if (i == 1) paste("RMSE <", t) else NULL) +  # Add titles only for the top row
        xlab(if(i==4) 'Time' else NULL) +
        ylab(if(j==1) out_name else NULL) +
        theme_minimal() 
      # Store the plot in the list
      plots[[length(plots) + 1]] <- p
    }
  }
  
  # Combine the plots into a grid
  final_plot <- wrap_plots(plots, ncol = length(thresholds))
  final_plot <- final_plot + 
    plot_annotation(
      title = main_title,
      subtitle = subtitle,
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5)
      )
    )
  
  print(final_plot)
}

search_space <- function(out_list, nout=4){
  all_results <- list()
  outs <- names(out_list)[1:nout]
  for (i in seq_along(outs)) {
    out_name <- outs[i]
    out <- out_list[[out_name]]
    dfx <- as_tibble(out$x_mat)
    colnames(dfx) <- c('p1', 'p2', 'rep')
    dfx$rmse <- out$y_rmse
    dfx$method <- out_name
    dfx$time <- 1:nrow(dfx)
    dfx$iter <- out$iter_idx
    all_results[[i]] <- dfx
  }
  return(do.call(rbind, all_results))
}

plot_threshold_search <- function(out_list, par_true, thresholds=c(15, 20, 25, 30)){
  outs <- names(out_list)[1:4]  # Assuming out_list is a named list
  config_info <- out_list$config
  main_title <- paste("Experiment", config_info$exp_id)
  subtitle <- paste(
    "Random seed =", config_info$exp_seed,
    "| Initial # parameters =", config_info$init_npar,
    "| Initial # reps =", config_info$nrep,
    "| Budget =", config_info$sim_budget,
    "| Grid size =", config_info$grid_npar,
    "| TS sample size =", config_info$nTS_samp
  )
  
  # Create a list to store all plots
  plots <- list()
  
  # Loop over each combination of thresholds and outs
  for (i in seq_along(outs)) {
    out_name <- outs[i]
    out <- out_list[[out_name]]
    dfx <- as_tibble(out$x_mat)
    colnames(dfx) <- c('p1', 'p2', 'rep')
    dfx$rmse <- out$y_rmse
    dfx$method <- out_name
    dfx$time <- 1:nrow(dfx)
    dfx$iter <- out$iter_idx
    for (j in seq_along(thresholds)) {
      t <- thresholds[j]
      dfx$within_threshold <- ifelse(dfx$rmse <= t, 1, 0)
      dfx$quality <- factor(dfx$within_threshold)
      # Generate the plot
      p <- ggplot(dfx, aes(x = p1, y = p2, color=quality))  +
        geom_point(size=3, alpha=.1) + 
        geom_point(data=data.frame(p1=par_true[1], p2=par_true[2], rep=par_true[3]), color='red', fill='red', size=2, shape=23) +
        scale_color_manual(values = c("grey", "green")) + 
        ggtitle(if (i == 1) paste("RMSE <", t) else NULL) +  # Add titles only for the top row
        xlab(if(i==4) 'Time' else NULL) +
        ylab(if(j==1) out_name else NULL) +
        theme_minimal() +
        theme(legend.position = "none") 
      # Store the plot in the list
      plots[[length(plots) + 1]] <- p
    }
  }
  
  # Combine the plots into a grid
  final_plot <- wrap_plots(plots, ncol = length(thresholds))
  
  print(final_plot)
}


plot_nbest_traj <- function(nbest, res_fixed, res_adaptive, yobs, nsteps = 100, legend_labels = c("Fixed", "Adaptive")){
  
  best_ix_fixed <- order(res_fixed$y_vec)[1:nbest]
  y_fixed_mat <- res_fixed$simout_mat[,best_ix_fixed]
  
  best_ix_adapt <-order(res_adaptive$y_vec)[1:nbest]
  y_adaptive_mat <- res_adaptive$simout_mat[,best_ix_adapt]
  
  
  df_traj <- as.data.frame(cbind(rep(0:nsteps, 2*nbest),
                                 c(as.numeric(y_fixed_mat), 
                                   as.numeric(y_adaptive_mat)),
                                 rep(1:(nbest*2), each = nsteps+1)))
  
  df_traj <- cbind(df_traj, factor(rep(legend_labels, each = (nsteps+1)*nbest)))
  df_true <- data.frame(tt = 0:nsteps, yobs = yobs)                         
  
  colnames(df_traj) <- c("time", "new_infections", "rep", "type")
  
  ggplot(data = df_traj) + 
    geom_line(mapping = aes(x = time, y = new_infections, color = type,
                            linetype = type, group = interaction(type, rep)), alpha = 0.5) + 
    geom_line(data = df_true, mapping = aes(x = tt, y = yobs), color = "black", linewidth = 1.5) +
    facet_wrap(~type) +
    theme_bw() +
    labs(title = paste("nbest = ", nbest))
}