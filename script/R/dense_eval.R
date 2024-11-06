## Dense evaluations for comparison

source("odin-sim.R")

nparam_dense <- 500
nrep <- 20

X_01 <- as.matrix(expand.grid(seq(0, 1, length.out = nparam_dense), 
                              seq(0, 1, length.out = nparam_dense)))
s <- 1:nrep
Xs_01 <- cbind(X_01[rep(1:nparam_dense^2, each = nrep), ], 
               rep(s, nparam_dense^2))

y_dense <- rep(NA, nparam_dense^2*nrep)
k <- 1
for (ii in 1:nparam_dense^2){
  for (seed in s){
    y_dense[k] <- run_sim_err(Xs_01[k, ])
    k <- k + 1
  }
}
d_dense <- data.frame(cbind(Xs_01, y_dense))
colnames(d_dense) <- c("beta", "gamma", "seed", "y")

save(d_dense, file = "../../data/dense_eval.RData")
