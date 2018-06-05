library(sl3)
library(origami)

Qbar0 <- function(A, W) {
  W1 <- W[, 1]
  W2 <- W[, 2]
  W3 <- W[, 3]
  W4 <- W[, 4]
  Qbar <- ifelse(W4 > 0, plogis(1 - W1^2 + 3 * W2 + A * (5 * W3^2 - 4.45)), plogis(-0.5 - W3 + 2 * W1 * W2 + A * (3 *
                                                                                                                    abs(W2) - 1.5)))
  return(Qbar)
}

g0 <- function(W) {
  W1 <- W[, 1]
  W2 <- W[, 2]
  W3 <- W[, 3]
  W4 <- W[, 4]

  plogis(0.25 * W1 - 0.1 * W2)
}

gen_data <- function(n = 1000, p = 4) {
  W <- matrix(rnorm(n * p), nrow = n)
  colnames(W) <- paste("W", seq_len(p), sep = "")
  A <- rbinom(n, 1, g0(W))
  u <- runif(n)
  Y <- as.numeric(u < Qbar0(A, W))
  Y0 <- as.numeric(u < Qbar0(0, W))
  Y1 <- as.numeric(u < Qbar0(1, W))
  d0 <- as.numeric(Qbar0(1, W) > Qbar0(0, W))
  Yd0 <- as.numeric(u < Qbar0(d0, W))
  data.frame(W, A, Y, Y0, Y1, Yd0, d0, blip = Qbar0(1, W) - Qbar0(0, W))
}

#Also contains conterfactuals, true rule and blip.
data_full <- gen_data(1000, 5)
data<-data_full[,1:7]

#Setup:
V=10
stratifyAY=TRUE
Wnodes=grep("^W", names(data), value = T)
Anode = "A"
Ynode = "Y"
Vnodes = Wnodes
Q_library = list("Lrnr_mean", "Lrnr_glm_fast")
g_library=list("Lrnr_mean", "Lrnr_glm_fast")
b_library=list("Lrnr_mean", "Lrnr_glm_fast", "Lrnr_glmnet","Lrnr_randomForest","Lrnr_xgboost")
gbounds=c(1e-4,1-1e-4)
Qbounds=c(1e-4,1-1e-4)
verbose = TRUE




