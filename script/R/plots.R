library(ggplot2)
library(GGally)

#' Title
#'
#' @param x_fixed_mat 
#' @param x_adaptive_mat 
#' @param save 
#'
#' @return
#' @export
#'
#' @examples
pairs_plot <- function(x_fixed_CRNGP_mat, 
                       x_adaptive_CRNGP_mat, 
                       x_fixed_hetGP_mat, 
                       x_adaptive_hetGP_mat,
                       save = F){
  
  df_fixed_CRNGP <- as.data.frame(x_fixed_CRNGP_mat)
  df_fixed_CRNGP$type <- "Fixed CRNGP grid"
  colnames(df_fixed_CRNGP)[1:3] <- c("beta", "gamma", "seed")
  
  df_adaptive_CRNGP <- as.data.frame(x_adaptive_CRNGP_mat)
  df_adaptive_CRNGP$type <- "Adaptive CRNGP grid"
  colnames(df_adaptive_CRNGP)[1:3] <- c("beta", "gamma", "seed")
  
  df_fixed_hetGP <- as.data.frame(x_fixed_hetGP_mat)
  df_fixed_hetGP$type <- "Fixed hetGP grid"
  colnames(df_fixed_hetGP)[1:3] <- c("beta", "gamma", "seed")
  
  df_adaptive_hetGP <- as.data.frame(x_adaptive_hetGP_mat)
  df_adaptive_hetGP$type <- "Adaptive hetGP grid"
  colnames(df_adaptive_hetGP)[1:3] <- c("beta", "gamma", "seed")
  
  df <- rbind(df_fixed_CRNGP, df_adaptive_CRNGP, df_fixed_hetGP, df_adaptive_hetGP)
  
  p <- ggpairs(df, 
               aes(color = type, fill = type, alpha = 0.5),
               columns = 1:2,
               diag = list(continuous = wrap("densityDiag", alpha = 0.5)),  # Density on diagonal
               # upper = list(continuous = wrap("density", alpha = 0.5)),  # 2D density on upper panels
               lower = list(continuous = wrap("density", alpha = 0.5))) +
               # diag = list(continuous = "barDiag")) +
    theme_minimal() +
    labs(title = "Evaluated points")
  
  if(save) ggsave("pairs.png")
  
  return(p)
  
}


#' Title
#'
#' @param nbest 
#' @param x_fixed_mat 
#' @param x_adaptive_mat 
#' @param y_fixed_vec 
#' @param y_adaptive_vec 
#' @param yobs 
#'
#' @return
#' @export
#'
#' @examples
plot_nbest_traj <- function(nbest, x_fixed_mat, x_adaptive_mat, y_fixed_vec, y_adaptive_vec, 
                            yobs, legend_labels = c("Fixed", "Adaptive")){
  
  x_fixed_best <- x_fixed_mat[order(y_fixed_vec)[1:nbest], ]
  x_adaptive_best <- x_adaptive_mat[order(y_adaptive_vec)[1:nbest], ]
  
  nsteps <- length(yobs) - 1
  
  # df_traj <- as.data.frame(matrix(NA, nrow = nsteps * nbest * 2), ncol = 3)
  
  y_fixed_mat <- y_adaptive_mat <- matrix(NA, ncol = nbest, nrow = nsteps+1)
  for(ii in 1:nbest){
    o_fixed <- run_sim(x_fixed_best[ii, ])
    o_adaptive <- run_sim(x_adaptive_best[ii, ])
    
    y_fixed_mat[, ii] <- o_fixed$n_SI
    y_adaptive_mat[, ii] <- o_adaptive$n_SI
  }
  
  df_traj <- as.data.frame(cbind(rep(0:nsteps, 2*nbest),
                                 c(as.numeric(y_fixed_mat), 
                                   as.numeric(y_adaptive_mat)),
                                 rep(1:(nbest*2), each = nsteps+1)))
  df_traj <- cbind(df_traj, factor(rep(legend_labels, each = (nsteps+1)*nbest)))
  df_true <- data.frame(tt = 0:nsteps, yobs = yobs)                         
  
  colnames(df_traj) <- c("time", "new_infections", "rep", "type")
  
  p <- ggplot(data = df_traj) + 
    geom_line(mapping = aes(x = time, y = new_infections, color = type,
                  linetype = type, group = interaction(type, rep)), alpha = 0.5) + 
    geom_line(data = df_true, mapping = aes(x = tt, y = yobs), color = "black", linewidth = 1.5) +
    facet_wrap(~type) +
    theme_bw() +
    labs(title = paste("nbest = ", nbest))
    
  
  return(p)
}


#' Title
#'
#' @param nbest 
#' @param y_fixed_native 
#' @param y_adaptive_native 
#' @param nerror_ticks 
#'
#' @return
#' @export
#'
#' @examples
plot_pct_traj <- function(nbest, y_fixed_native, y_adaptive_native, 
                          nerror_ticks = 20){
  
  y_df <- data.frame(y = c(y_fixed_native[order(y_fixed_native)][1:nbest], 
                           y_adaptive_native[order(y_adaptive_native)][1:nbest]),
                   type = rep(c("Fixed", "Adaptive"), each = nbest))

  error_ticks <- seq(0, max(y_df$y), length.out = nerror_ticks)
  pctraj_adapt <- pctraj_random <- rep(NA, nerror_ticks)
  
  for (ii in 1:nerror_ticks){
    pctraj_adapt[ii] <- mean(y_df$y[y_df$type == "Adaptive"] < error_ticks[ii])
    pctraj_random[ii] <- mean(y_df$y[y_df$type == "Fixed"] < error_ticks[ii])
  }
  
  plot(error_ticks, pctraj_adapt, pch = 19, col = alpha("#F8766D", 01), type = "b", 
       xlab = "sum-squared-error", ylab = "% best traj", main = paste("nbest = ", nbest))
  points(error_ticks, pctraj_random, pch = 19, col = alpha("#00BA38", 1), type = "b")
  legend("topleft", c("Adaptive", "Fixed"), pch = 19, lty = 1, col = c("#F8766D", "#00BA38"))
}

#' Title
#'
#' @param nbest 
#' @param y_fixed_native 
#' @param y_adaptive_native 
#' @param nerror_ticks 
#'
#' @return
#' @export
#'
#' @examples
plot_pct_traj2 <- function(nbest, y1, y2, labs, 
                          nerror_ticks = 20){
  
  y_df <- data.frame(y = c(y1[order(y1)][1:nbest], 
                           y2[order(y2)][1:nbest]),
                     type = rep(labs, each = nbest))
  
  error_ticks <- seq(0, max(y_df$y), length.out = nerror_ticks)
  pctraj_1 <- pctraj_2 <- rep(NA, nerror_ticks)
  
  for (ii in 1:nerror_ticks){
    pctraj_1[ii] <- mean(y_df$y[y_df$type == labs[1]] < error_ticks[ii])
    pctraj_2[ii] <- mean(y_df$y[y_df$type == labs[2]] < error_ticks[ii])
  }
  
  plot(error_ticks, pctraj_1, pch = 19, col = alpha("#F8766D", 01), type = "b", 
       xlab = "sum-squared-error", ylab = "% best traj", main = paste("nbest = ", nbest))
  points(error_ticks, pctraj_2, pch = 19, col = alpha("#00BA38", 1), type = "b")
  legend("topleft", labs, pch = 19, lty = 1, col = c("#F8766D", "#00BA38"))
}


#' Title
#'
#' @param nbest 
#' @param y_fixed_native 
#' @param y_adaptive_native 
#'
#' @return
#' @export
#'
#' @examples
err_hist <- function(nbest, y_fixed_native, y_adaptive_native){
  
  y_df <- data.frame(y = c(y_fixed_native[order(y_fixed_native)][1:nbest], 
                           y_adaptive_native[order(y_adaptive_native)][1:nbest]),
                     type = rep(c("Fixed", "Adaptive"), each = nbest))
  
  p <- ggplot(y_df, aes(x = y, fill = type)) + 
              geom_histogram(alpha = 0.5) + 
              # geom_vline(aes(xintercept = 0)) +
              xlim(min(y_df$y)-1e4, max(y_df$y)+1e4) +
              facet_wrap(~type) + 
              xlab("sum sqaured error") +
              theme_bw() +
    labs(title = paste("nbest = ", nbest))
  
  return(p)
}

#' Title
#'
#' @param nbest 
#' @param rep_id 
#' @param x_fixed_mat 
#' @param x_adaptive_mat 
#' @param y_fixed_vec 
#' @param y_adaptive_vec 
#' @param d_dense 
#'
#' @return
#' @export
#'
#' @examples
contour_compare <- function(nbest, rep_id, x_fixed_mat, x_adaptive_mat, y_fixed_vec, y_adaptive_vec, d_dense){
  
  x_fixed_best <- x_fixed_mat[order(y_fixed_vec)[1:nbest], ]
  x_fixed_best <- x_fixed_best[x_fixed_best[, 3] == rep_id, ]
  if(!is.matrix(x_fixed_best)) x_fixed_best <- matrix(x_fixed_best, nrow = 1)
  
  x_adaptive_best <- x_adaptive_mat[order(y_adaptive_vec)[1:nbest], ]
  x_adaptive_best <- x_adaptive_best[x_adaptive_best[, 3] == rep_id, ]
  if(!is.matrix(x_adaptive_best)) x_adaptive_best <- matrix(x_adaptive_best, nrow = 1)
  
  ggplot(d_dense) +
    geom_contour(aes(x = beta, y = gamma, z = y, color = ..level..), bins = 40) +  
    geom_point(data = data.frame(x_adaptive_best[, 1:2]), 
               mapping = aes(x = X1, y = X2), color = "blue", size = 3) +
    
    geom_point(data = data.frame(x_fixed_best[, 1:2]), 
               mapping = aes(x = X1, y = X2), color = "red", size = 3) +
    
    # scale_color_brewer(palette = "Greens") +
    ggtitle("Blue = adaptive, Red = random")
  theme_minimal()
}


## TODO: area under the curve

area_under_curve <- function(nbests, y1, y2, labs,
                             nerror_ticks = 20){
  
  area_vec <- rep(NA, length(nbests))
  
  l <- 1
  for (nbest in nbests){
    y_df <- data.frame(y = c(y1[order(y1)][1:nbest], 
                             y2[order(y2)][1:nbest]),
                       type = rep(labs, each = nbest))
    
    error_ticks <- seq(0, max(y_df$y), length.out = nerror_ticks)
    pctraj_1 <- pctraj_2 <- rep(NA, nerror_ticks)
    
    for (ii in 1:nerror_ticks){
      pctraj_1[ii] <- mean(y_df$y[y_df$type == labs[1]] < error_ticks[ii])
      pctraj_2[ii] <- mean(y_df$y[y_df$type == labs[2]] < error_ticks[ii])
    }
  
    area_vec[l] <- sum(pctraj_1 - pctraj_2)
    l <- l + 1
  }
  
  barplot(area_vec, names.arg = nbests, 
          ylab = "area under curve", xlab = "nbest trajectories",
          main = paste(labs[1], " - ", labs[2]))
  
}
