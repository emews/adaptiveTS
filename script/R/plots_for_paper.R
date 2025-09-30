library(hetGP)
library(scales)
library(data.table)
library(lhs)
library(ggplot2)
library(mvtnorm)
library(tidyr)
library(dplyr)
library(plotly)
library(patchwork)

source("odin-sim.R")
source("adaptive_TS.R")
source("adaptive_TS_nonCRN.R")
source("plots.R")
source("analysis_helper_abby.R")

doSave <- TRUE

## ggthemes for 3 panel plots
theme_3 <- theme(
    strip.text = element_text(size = 12, face = "bold"),         # Facet title
    axis.title = element_text(size = 12),                        # Axis titles
    axis.text = element_text(size = 12),                         # Axis tick labels
    legend.title = element_text(size = 12),                      # Legend title
    legend.text = element_text(size = 12),                       # Legend labels
    strip.background = element_rect(fill = "white", color = "black")
  )
  

## ============================================================
## Experiment meta data and true trajectory

beta_true  <- 0.7
gamma_true <- 0.2
seed_true <- 50

# exp_path  <- "../../emews/experiments/all_methods_sweep/results"
exp_path <- "../../emews/experiments/exp_gt_beta0.7_gamma0.2_seed50/results/"
plot_path <- paste0("../../plots/experiments/exp_gt_beta", beta_true, "_gamma", gamma_true, "_seed50/")
config <- fread("../../emews/data/upfs/config_dedup.csv")
config[, exp_id := as.numeric(gsub("exp_", "", exp_id))]
par_true <- c(beta_true, gamma_true, seed_true)

## ============================================================
# Get true trajectory 
xtrue <- par_true[1:2]
strue <- par_true[3]
ytrue_sim <- run_sim(c(xtrue, strue))
ytrue <-  ytrue_sim$n_SI

sim_seeds <- 1:50
ytrue_sims <- list()
ytrues <- list()
i <- 1
for (s in sim_seeds){
  sim <- run_sim(c(xtrue, s))
  sim$seed <- s
  ytrue_sims[[i]] <- sim
  ytrues[[i]] <- sim$n_SI
  i <- i + 1
}

gt_sims <- do.call(rbind, ytrue_sims)
gt_sims$seed <- as.factor(gt_sims$seed)

# First, compute the mean_SI by step
gt_summary <- gt_sims %>%
  group_by(step) %>%
  summarise(mean_SI = mean(n_SI), .groups = "drop")

# Merge the summary back into gt_sims by 'step'
gt_agg <- merge(gt_sims, gt_summary, by = "step")


## look at replicates at ground thuth parameter
ggplot(gt_sims) + 
  geom_line(aes(x=step, y=n_SI, group=seed), color='grey', size=.4) + 
  geom_line(data=gt_sims %>% filter(seed == par_true[3]), aes(x=step, y=n_SI, group=seed), color='forestgreen', size=1) +
  #geom_line(data=gt_sims %>% filter(seed %in% unique(config$sim_seed)), aes(x=step, y=n_SI, color=seed), size=1)+
  geom_line(data=gt_sims %>% group_by(step) %>% summarise(n_SI = mean(n_SI)), aes(x=step, y=n_SI), color='black', linewidth=.75) +
  xlab('Time') + 
  ylab('New Infections') +
  theme_minimal()

## read all experiments
# all_out_lists <- lapply(1:max(config$exp_id), load5outputs, exp_path=exp_path, config=config)
# saveRDS(all_out_lists, file = paste0(exp_path, "combined_out_lists.rds"))
all_out_lists <- readRDS(paste0(exp_path, "combined_out_lists.rds"))

## ============================================================
## Plot trajectories colored by RMSE bins for aCRNGP method
## ============================================================

out_list <- all_out_lists[[10]]
out <- out_list$adaptive_crn
yobs <- ytrue_sim$n_SI

df_true <- data.frame(time= (1:length(yobs)) - 1, yobs = yobs)
df_traj <- as.data.frame(t(out$simout_mat))
colnames(df_traj) <- 1:ncol(df_traj)
df_traj$y_rmse <- out$y_rmse
df_traj$traj <- 1:nrow(df_traj)

#df_traj$time <- (1:nrow(df_traj)) - 1
df_traj <- pivot_longer(df_traj, 
               cols = -c(y_rmse, traj), 
               names_to = "time", 
               values_to = "new_infections")
df_traj$time <- as.numeric(df_traj$time)

label_rmse_bin <- function(y_vals, thresholds) {
  # Make sure thresholds is sorted ascending
  thresholds <- sort(thresholds)
  
  # For each y, pick the first t in thresholds such that y < t.
  # If none exist, assign "RMSE >= max(thresholds)".
  
  sapply(y_vals, function(y) {
    larger <- thresholds[y < thresholds]
    if (length(larger) > 0) {
      paste0("RMSE < ", larger[1])
    } else {
      paste0("RMSE >= ", thresholds[length(thresholds)])
    }
  })
}

df_traj$threshold <- label_rmse_bin(df_traj$y_rmse, c(15, 20, 25, 30))

my_colors <- c(
  "RMSE < 15" = "#009E73", 
  "RMSE < 20" = "#F0E442",  
  "RMSE < 25" = "#56B4E9",
  "RMSE < 30" = "#D55E00"
)

p_gt <- 
  ggplot() +
    geom_line(
      data = df_traj %>% filter(threshold == 'RMSE >= 30'), 
      aes(x = time, 
          y = new_infections, 
          group = traj
      ),
      alpha = 0.5,
      color='grey'
    ) +
    geom_line(
      data = df_traj %>% filter(threshold == 'RMSE < 30'), 
      aes(x = time, 
          y = new_infections, 
          group = traj, 
          color = factor(threshold)   # ensure it's treated as discrete
      ),
    ) +
    geom_line(
      data = df_traj %>% filter(threshold == 'RMSE < 25'), 
      aes(x = time, 
          y = new_infections, 
          group = traj, 
          color = factor(threshold)   # ensure it's treated as discrete
      ),
    ) +
    geom_line(
      data = df_traj %>% filter(threshold == 'RMSE < 20'), 
      aes(x = time, 
          y = new_infections, 
          group = traj, 
          color = factor(threshold)   # ensure it's treated as discrete
      ),
    ) +
    geom_line(
      data = df_traj %>% filter(threshold == 'RMSE < 15'), 
      aes(x = time, 
          y = new_infections, 
          group = traj, 
          color = factor(threshold)   # ensure it's treated as discrete
      ),
    ) +
    geom_point(
      data    = df_true, 
      mapping = aes(x = time, y = yobs),
      color = "black", 
      size = 1
    ) +
    scale_color_manual(
      name   = "Threshold", 
      values = my_colors
    ) +
    theme_bw() +
    labs(
      x = "Time",
      y = "New infections"
    ) +
    theme(
      legend.position = "inside",
      legend.position.inside = c(0.8, 0.8),
      legend.background = element_blank()
    )
p_gt
if(doSave){
  ggsave(paste0(plot_path, "trajectories_found_aCRN_example_", beta_true, "_", gamma_true, ".png"), p_gt, width=6, height=4)
}

## ============================================================
## Number of trajectories sampled by each method
## ============================================================

# count the number of trajectories found under each threshold for each experiment
df_threshold <- do.call(
  rbind,
  lapply(c(300, 500, 700), function(b) {
    do.call(
      rbind,
      lapply(all_out_lists, function(x) {
        threshold_counts(x, sim_budget = b)
      })
    )
  })
)

df_threshold <- merge(df_threshold, config[,sim_budget:=NULL])
method_names <- c(
  "fixed_crn" = "fCRN",      
  "adaptive_crn"  = "aCRN",
  "fixed_het"      = "fHet",
  "adaptive_het"   = "aHet",
  "fixed_crn_seed" = "fgCRN"
)

dt_threshold <- as.data.table(df_threshold)
dt_threshold[, method := method_names[method]]
dt_threshold$method <- factor(dt_threshold$method, levels=c("fHet", "aHet", "aCRN", "fgCRN", "fCRN"))

## overall comparison of methods across thresholds and budgets and experiments
dt_threshold$pct_traj <- df_threshold$num_trajectories/df_threshold$sim_budget
dt_threshold$sim_budget_lab <- paste0("N[max]==", dt_threshold$sim_budget)

dt_exp <- dt_threshold %>% 
            group_by(threshold, sim_budget, sim_budget_lab, method, init_npar, nrep, grid_npar, nTS_samp, exp_id) %>%
            summarise(mean_traj=mean(num_trajectories))
dt_exp$pct_traj <- dt_exp$mean_traj/dt_exp$sim_budget
dt_exp$method <- factor(dt_exp$method, levels=c("aCRN", "fHet", "aHet", "fgCRN", "fCRN"))

# overall boxplot of proportion of budget used to find trajectories under each threshold
p_pct_traj <- ggplot(dt_exp, aes(x=factor(threshold), y=pct_traj, fill=method)) + 
  #geom_bar(stat="summary", position="dodge") + 
  #geom_errorbar(stat='summary', position="dodge") +
                    geom_boxplot() +
                    xlab('Threshold (RMSE)') +
                    ylab('Prop. of budget under threshold') +
                    facet_wrap(vars(sim_budget_lab), labeller = label_parsed, ncol = 3) +
                    # facet_wrap(~sim_budget, ncol=3, labeller = as_labeller(sim_budget = function(x) paste("N[max]~`=`~", x)))+
                    theme_bw() +
                    theme_3
p_pct_traj
if(doSave){
  ggsave(paste0(plot_path, "pct_traj_boxplot_", beta_true, "_", gamma_true, ".png"), p_pct_traj, width=12, height=4)
}

# Aggregate over designs 
hyper_designs <- df_threshold %>% 
  select(sim_budget, init_npar, nrep, grid_npar, nTS_samp) %>%
  distinct() %>%
  mutate(
    exp_design = row_number(),
    initial_sims = init_npar * nrep,
    sampling_size = grid_npar * nTS_samp)

# find the best performing design for each method, threshold, budget
df_threshold_best <- dt_exp %>%
  group_by(threshold, sim_budget, method) %>%
  slice_max(mean_traj, with_ties = FALSE) %>%
  ungroup()
# Now get ALL seeds/experiments that match these best design combinations
df_threshold_best_reps <- dt_threshold %>%
  # Join with the best designs to filter only those combinations
  inner_join(
    df_threshold_best %>% select(threshold, sim_budget, method, init_npar, nrep, grid_npar, nTS_samp),
    by = c("threshold", "sim_budget", "method", "init_npar", "nrep", "grid_npar", "nTS_samp")
  ) %>%
  # Calculate pct_traj for each individual experiment/seed
  mutate(pct_traj = num_trajectories / sim_budget)

# summarize over experiments for each method, threshold, budget
df_summary <- merge(df_threshold_best_reps, hyper_designs) %>%
  group_by(sim_budget, threshold, method, exp_design, sim_budget_lab) %>%
  # summarise(pct_traj=mean(pct_traj)) %>%
  # group_by(sim_budget, threshold, method) %>%
  summarize(
    median_prop = median(pct_traj),
    iqr_lower   = quantile(pct_traj, 0.25),
    iqr_upper   = quantile(pct_traj, 0.75),
    mean_prop = mean(pct_traj),
    se_lower = mean(pct_traj) - sd(pct_traj)/sqrt(n()),
    se_upper = mean(pct_traj) + sd(pct_traj)/sqrt(n())
  )
df_summary$method <- factor(df_summary$method, levels=c("aCRN", "fHet", "aHet", "fgCRN", "fCRN"))

p_pct_traj_best <- ggplot(df_summary, aes(x = threshold, y = median_prop, color = method)) +
                          geom_line(size = 1) +
                          geom_point(size = 2) +
                          geom_ribbon(aes(ymin = iqr_lower, ymax = iqr_upper, fill = method),
                                      alpha = 0.2, color = NA) +
                          facet_wrap(vars(sim_budget_lab), labeller = label_parsed, ncol = 3) +
                          scale_x_continuous(breaks = c(15, 20, 25, 30)) +
                          labs(
                            x = "Threshold (RMSE)",
                            y = "Median Proportion of Budget Under Threshold",
                            color = "Method",
                            fill = "Method"
                          ) +
                          theme_bw() + theme_3

p_pct_traj_best
if(doSave){
  ggsave(paste0(plot_path, "pct_traj_bestdesign_line_", beta_true, "_", gamma_true, ".png"), p_pct_traj_best, width=12, height=4)
}

# boxplots showing effect of size of initial design
df_pairs <- dt_threshold %>% filter(threshold==30)
df_pairs$initial_sims <- df_pairs$init_npar * df_pairs$nrep
df_pairs$init_cost_lab  <- paste0("N[0]==", df_pairs$initial_sims)
df_pairs$method <- factor(df_pairs$method, levels=c("aCRN", "fHet", "aHet", "fgCRN", "fCRN"))

p_init_sims <- ggplot(df_pairs) + 
  geom_boxplot(aes(x = method, y = pct_traj, fill = factor(initial_sims, levels = c(200, 100, 50)))) +
  facet_wrap(~sim_budget_lab, labeller = label_parsed) +
  labs(
    x = "Method",
    y = "Proportion under threshold RMSE < 30",
    fill = expression(N[0]),
    color = expression(N[0])
  ) +
  theme_bw() + 
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12)
  ) + theme_3 +
  scale_fill_manual(
    values = c("50" = "#2a68a7ff", "100" = "#6487a5ff", "200" = "#56B1F7")
  ) 
  # scale_color_manual(
  #   values = c("50" = "#2a68a7ff", "100" = "#6487a5ff", "200" = "#56B1F7")
  # )
p_init_sims
if(doSave){
  ggsave(paste0(plot_path, "pct_traj_boxplot_initdesign_", beta_true, "_", gamma_true, ".png"), p_init_sims, width=10, height=4)
}


## ============================================================
## Time to solution
## ============================================================

# for each method, each threshold, each experiment, compute AUC
df_auc <- do.call(
  rbind,
  lapply(all_out_lists, function(x) {
         threshold_budget_analysis(x, nout=5)$agg_auc
  })
)
df_auc <- df_auc %>% 
  mutate(auc_scaled=auc/(sim_budget**2))

# for each method, each threshold, each experiment, compute number of trajectories
# as a function of continuous budget
df_budget <- do.call(
  rbind,
  lapply(all_out_lists, function(x) {
         threshold_budget_analysis(x, nout=5)$budget_df
  })
)

# best design for each method, threshold, budget
df_best <- dt_exp %>%
  group_by(threshold, sim_budget, method) %>%
  slice_max(mean_traj, with_ties = FALSE) %>%
  ungroup()


df_auc <- merge(df_auc, config)
df_auc$method <- method_names[df_auc$method]
df_budget$method <- method_names[df_budget$method]

# for each method, threshold, budget = 700, get best design according to 
# number of trajectories found on average across exp_seeds and the collect
# all seeds AUC for that design

# Explicitly remove exp_id from one of the dataframes before merging
df_auc_best <- merge(
  df_auc,
  df_best %>% select(-exp_id),
  by = c("threshold", "sim_budget", "method", "init_npar", "nrep", "grid_npar", "nTS_samp")
)

# ggplot(df_auc_best, aes(x=factor(threshold), y=auc, color=method, fill=method)) + 
#   geom_boxplot() + 
#   facet_wrap(~sim_budget, labeller = labeller(sim_budget = function(x) paste('budget=', x)), scale='free_y')+
#   theme_minimal()

df_auc_best_agg <- df_auc_best %>% 
                group_by(threshold, sim_budget, method) %>% 
                summarise(mean_auc=mean(auc), rauc=mean(auc_scaled))

# ggplot(df_auc_best_agg, aes(x=threshold, y=mean_auc, color=method, group=method)) + 
#   geom_point()+
#   geom_line() +
#   theme_minimal()

df_budget_best <- merge(df_budget, df_auc_best, by=c('method', 'threshold', 'exp_id'))
df_budget_best <- df_budget_best[order(df_budget_best$time),]
df_budget_best_agg <- df_budget_best %>% 
                          group_by(threshold, method, time) %>% 
                          summarise(mean_found=mean(trajectories_found))

# ggplot(df_budget_best_agg %>% filter(threshold==20, time<=700)) + 
#   geom_line(aes(x=time, y=mean_found, color=method, group=method)) + 
#   theme_minimal()

df_auc_best_rep <- df_budget_best %>% filter(time %in% c(300, 500, 700))
df_auc_best_rep$sim_budget_lab <- paste0("N[max]==", df_auc_best_rep$time)
df_auc_best_rep$method <- factor(df_auc_best_rep$method, levels=c("aCRN", "fHet", "aHet", "fgCRN", "fCRN"))

pt1 <- ggplot(df_auc_best_rep, aes(x=factor(threshold), y=auc_scaled, fill=method)) + 
  geom_boxplot() + 
  facet_wrap(vars(sim_budget_lab), labeller = label_parsed, ncol = 3) +
  theme_bw() +
  labs(
    x = "Threshold (RMSE)",
    y = "rAUC",
    color = "Method",
    fill = "Method"
  ) +
  theme_3 
pt1
if(doSave){
  ggsave(paste0(plot_path, "auc_boxplot_bestdesign_", beta_true, "_", gamma_true, ".png"), pt1, width=12, height=4)
}



## only for the main paper experiment

t<-15
exp_ids <- c(90,110) # 310,35, exp_id
label_data = df_auc %>% filter(method == "aCRN", threshold==t, exp_id %in% exp_ids)
label_data$time = 700
label_data$text_pos = c(225, 250)

plot_data <- df_budget %>%
  filter(exp_id %in% exp_ids, threshold==15, method == "aCRN")
plot_data$method <- factor(plot_data$method, levels=c("aCRN", "fHet", "aHet", "fgCRN", "fCRN"))


# Plot with annotation
fig <- ggplot(plot_data) + 
  geom_line(aes(x = time, y = trajectories_found, color = factor(exp_id)), size=.75) +
  geom_text(
    data = label_data,
    aes(x = time, y = text_pos, 
        #label =paste0("rAUC=", sprintf("%.3f", auc_scaled)), 
        label=sprintf("%.3f", auc_scaled)),
        # color = factor(exp_id),
    hjust = -0.1,
    fontface="bold",
    show.legend=FALSE
  ) +
  theme_bw() +
  xlim(min(plot_data$time), max(plot_data$time) + 40)+
  labs(
    x = "Budget",
    y = paste0("Trajectories found with RMSE < ",t),
    color = "Experiment"
  ) +
  scale_color_manual(
    name = "Experiment",
    values = c("90" = "salmon", "110" = "cyan3"),  # Keep original colors
    labels = c("90" = "Configuration 1", "110" = "Configuration 2")
  ) +
 theme(
    legend.position = c(0.15, 0.95),
    legend.justification = c("left", "top"),
  legend.background = element_rect(fill = "white", color = NA),
    legend.box.background = element_blank()
  )
fig
if(doSave){
  ggsave(paste0(plot_path, "auc_pathway_example_", beta_true, "_", gamma_true, ".png"), fig, width=6, height=4)
}