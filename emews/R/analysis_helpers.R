library(data.table)
library(dplyr)
library(tidyr)
library(readr)

combine_results <- function(out_file, 
                            experiment_name='experiment', 
                            summary_file=NULL){
  out <- readRDS(out_file)
  if (length(out) ==2){
    out <- out$out
  }
  x_mat <- do.call(rbind, out$X_list)
  y_vec <- unlist(out$y_list)
  y_native <- unlist(out$ynative_list)
  iter_idx=unlist(mapply(rep, seq_along(out$X_list), sapply(out$X_list, nrow)))
  
  out_df <- as_tibble(x_mat)
  colnames(out_df) <- c('s2x', 'sah', 'ba', 'seed')
  out_df$y_std <- y_vec
  out_df$y <- y_native
  out_df$iter <- iter_idx
  out_df$budget <- 1:nrow(out_df)
  out_df$experiment <- experiment_name
  
  if(!is.null(summary_file)){
    summary <- read.csv(summary_file)
    # some previous analysis indicates that instance order corresponds to budget order
    # also rename something
    summary <- summary %>% arrange(instance) %>% rename(summary_seed=seed)
    out_df <- cbind(out_df, summary)
  }
  return(out_df)
}


prepare_trajectories <- function(output_file, 
                                 gt_h_file = '../data/calibration_targets/dt.chicago.hosp.csv',
                                 gt_d_file = '../data/calibration_targets/dt.chicago.deaths.csv',
                                 start_date = as.Date('2020-03-17'),
                                 end_date = as.Date("2020-06-13"),
                                 melted = TRUE) {
  
  # Load simulation data
  sim_df <- read_csv(output_file, col_types = cols()) %>%
    select(tick, hosp_r_count, hosp_icu_r_count, hosp_d_count, hosp_icu_d_count, icu_r_count, icu_d_count, dead_count) %>%
    mutate(date = as.Date("2020-03-16") + (tick / 24),
           total_hosp_sim = hosp_r_count + hosp_icu_r_count + hosp_d_count + hosp_icu_d_count + icu_r_count + icu_d_count)

# Load and filter ground truth data
gt_d_df <- read_csv(gt_d_file, col_types = cols()) %>%
  filter(date >= start_date & date <= end_date) %>%
  drop_na(deaths)

gt_h_df <- read_csv(gt_h_file, col_types = cols()) %>%
  filter(date >= start_date & date <= end_date) %>%
  drop_na(tot.hosp)

# Merge simulation and ground truth
df_all <- left_join(gt_h_df, sim_df, by = "date") %>%
  left_join(gt_d_df, by = "date") %>%
  mutate(err = abs(tot.hosp - total_hosp_sim) / tot.hosp + abs(deaths - dead_count) / deaths)

if (!melted) {
  return(df_all)
}

# Reshape into long format for plotting
df_melted <- df_all %>%
  pivot_longer(cols = c(total_hosp_sim, tot.hosp, dead_count, deaths),
               names_to = "type_metric",
               values_to = "val") %>%
  mutate(metric = ifelse(type_metric %in% c("tot.hosp", "total_hosp_sim"), "hosp", "deaths"),
         type = ifelse(type_metric %in% c("tot.hosp", "deaths"), "ground_truth", "simulation")) %>%
  select(date, tick, metric, type, val)  # Drop unnecessary columns

return(df_melted)
}
