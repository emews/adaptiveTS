library(data.table)

# output_file <- "/lcrc/project/EMEWS/improv/ncollier/repos/snl_brave_citycovid/experiments/stein_mcmc_01212025_1.0/instances/instance_0/output/counts_r0_0.csv"
# gt_h_file <- "/lcrc/project/EMEWS/afadikar/git/CityCOVID_DA/data/dt.chicago.hosp.csv"
# gt_d_file <- "/lcrc/project/EMEWS/afadikar/git/CityCOVID_DA/data/dt.chicago.deaths.csv"

#' Objective on hospitalization
#'
#' @param output_file path to the CityCOVID count file
#' @param gt_h_file path to the chicago ground truth hospitalization
#' @param start_date beginning date of calibration period
#' @param end_date end date of calibration period
#'
#' @return sum of squared difference between simulated and observed trajectory

obj_h <- function(output_file, gt_h_file, 
                  start_date = as.Date('2020-03-17'), 
                  end_date = as.Date("2020-06-13")){
  
  ## process simulation
  sim_df <- data.table::fread(output_file,
                              select = c("tick", "hosp_r_count", "hosp_icu_r_count",
                                         "hosp_d_count", "hosp_icu_d_count",
                                         "icu_r_count", "icu_d_count"))
  sim_df[, date := as.IDate("2020-03-16") + (tick / 24)]
  sim_df[, total_hosp_sim := rowSums(.SD), .SDcols = patterns("count$")]
  
  ## process ground truth
  gt_df <- data.table::fread(gt_h_file)
  gt_df <- gt_df[(date >= start_date) & (date <= end_date)]
  gt_df <- gt_df[!is.na(tot.hosp)]
  
  ## merge
  df_all <- merge(sim_df, gt_df, by = "date", all.y = TRUE)
  
  return(sum((df_all$tot.hosp - df_all$total_hosp_sim)^2))
  
}

#' Objective on deaths
#'
#' @param output_file path to the CityCOVID count file
#' @param gt_d_file path to the chicago ground truth deaths
#' @param start_date beginning date of calibration period
#' @param end_date end date of calibration period
#'
#' @return sum of squared difference between simulated and observed trajectory

obj_d <- function(output_file, gt_d_file, 
                  start_date = as.Date('2020-03-17'), 
                  end_date = as.Date("2020-06-13")){
  
  ## process simulation
  sim_df <- data.table::fread(output_file,
                              select = c("tick", "dead_count"))
  sim_df[, date := as.IDate("2020-03-16") + (tick / 24)]
  
  ## process ground truth
  gt_df <- data.table::fread(gt_d_file)
  gt_df <- gt_df[(date >= start_date) & (date <= end_date)]
  gt_df <- gt_df[!is.na(deaths)]
  
  ## merge
  df_all <- merge(sim_df, gt_df, by = "date", all.y = TRUE)
  
  return(sum((df_all$dead_count - df_all$deaths)^2))
  
}



#' Objective on hospitalizations and deaths
#'
#' @param output_file path to the CityCOVID count file
#' @param gt_h_file path to the chicago ground truth hospitalization
#' @param gt_d_file path to the chicago ground truth deaths
#' @param start_date beginning date of calibration period
#' @param end_date end date of calibration period
#'
#' @return sum of relative errors between simulated and observed trajectories

obj_dh <- function(output_file, gt_h_file, gt_d_file, 
                   start_date = as.Date('2020-03-17'), 
                   end_date = as.Date("2020-06-13")){
  
  ## process simulation
  sim_df <- data.table::fread(output_file,
                              select = c("tick", "hosp_r_count", "hosp_icu_r_count",
                                         "hosp_d_count", "hosp_icu_d_count",
                                         "icu_r_count", "icu_d_count", "dead_count"))
  sim_df[, date := as.IDate("2020-03-16") + (tick / 24)]
  sim_df[, total_hosp_sim := hosp_r_count + hosp_icu_r_count +
           hosp_d_count + hosp_icu_d_count + icu_r_count + icu_d_count]
  
  ## process ground truth
  gt_d_df <- data.table::fread(gt_d_file)
  gt_h_df <- data.table::fread(gt_h_file)
  
  gt_d_df <- gt_d_df[(date >= start_date) & (date <= end_date)]
  gt_h_df <- gt_h_df[(date >= start_date) & (date <= end_date)]
  
  gt_d_df <- gt_d_df[!is.na(deaths)]
  gt_h_df <- gt_h_df[!is.na(tot.hosp)]
  
  ## merge
  df_all <- merge(sim_df, gt_h_df, by = "date", all.y = TRUE)
  df_all <- merge(df_all, gt_d_df, by = "date", all.x = TRUE)
  
  df_all[, err:= abs(tot.hosp - total_hosp_sim)/tot.hosp + abs(deaths - dead_count)/deaths]
  
  return(sum(df_all$err))
  
}
