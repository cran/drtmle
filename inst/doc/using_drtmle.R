## ---- echo = FALSE, message = FALSE--------------------------------------
library(drtmle)
set.seed(1234)
N <- 1e6
W <- data.frame(W1 = runif(N), W2 = rbinom(N, 1, 0.5))
Y1 <- rbinom(N, 1, plogis(W$W1 + W$W2))
Y0 <- rbinom(N, 1, plogis(W$W1))
EY1 <- mean(Y1)
EY0 <- mean(Y0)

## ------------------------------------------------------------------------
set.seed(1234)
n <- 200
W <- data.frame(W1 = runif(n), W2 = rbinom(n, 1, 0.5))
A <- rbinom(n, 1, plogis(W$W1 + W$W2))
Y <- rbinom(n, 1, plogis(W$W1 + W$W2*A))

## ------------------------------------------------------------------------
glm_fit_uni <- drtmle(W = W, A = A, Y = Y, family = binomial(),
                      glm_g = "W1 + W2", glm_Q = "W1 + W2*A",
                      glm_gr = "Qn", glm_Qr = "gn", stratify = FALSE)
glm_fit_uni

glm_fit_biv <- drtmle(W = W, A = A, Y = Y, family = binomial(),
                      glm_g = "W1 + W2", glm_Q = "W1 + W2*A",
                      glm_gr = "Qn", glm_Qr = "gn", stratify = FALSE,
                      reduction = "bivariate")
glm_fit_biv

## ---- echo = FALSE, message = FALSE--------------------------------------
library(SuperLearner)

## ---- message = FALSE, cache = TRUE--------------------------------------
sl_fit <- drtmle(W = W, A = A, Y = Y, family = binomial(),
                 SL_g = c("SL.glm", "SL.step.interaction","SL.mean"),
                 SL_Q = c("SL.glm", "SL.step.interaction","SL.mean"),
                 SL_gr = c("SL.glm", "SL.earth", "SL.mean"),
                 SL_Qr = c("SL.glm", "SL.earth", "SL.mean"),
                 stratify = FALSE)
sl_fit

## ------------------------------------------------------------------------
SL.npreg

## ---- cache = TRUE, message = FALSE--------------------------------------
npreg_fit <- drtmle(W = W, A = A, Y = Y, family = binomial(),
                    SL_g = "SL.npreg", SL_Q = "SL.npreg",
                    SL_gr = "SL.npreg", SL_Qr = "SL.npreg",
                    stratify = TRUE)
npreg_fit

## ------------------------------------------------------------------------
ci(npreg_fit)

## ------------------------------------------------------------------------
ci(npreg_fit, contrast = c(1,-1))

## ------------------------------------------------------------------------
riskRatio <- list(f = function(eff){ log(eff) },
                  f_inv = function(eff){ exp(eff) },
                  h = function(est){ est[1]/est[2] },
                  fh_grad =  function(est){ c(1/est[1],-1/est[2]) })
ci(npreg_fit, contrast = riskRatio)

## ------------------------------------------------------------------------
logitMean <- list(f = function(eff){ qlogis(eff) },
                  f_inv = function(eff){ plogis(eff) },
                  h = function(est){ est[1] },
                  fh_grad = function(est){ c(1/(est[1] - est[1]^2), 0) })
ci(npreg_fit, contrast = logitMean)

## ------------------------------------------------------------------------
wald_test(npreg_fit, null = c(0.5, 0.6))

## ------------------------------------------------------------------------
wald_test(npreg_fit, contrast = c(1, -1))
wald_test(npreg_fit, contrast = c(1, -1), null = 0.05)

## ------------------------------------------------------------------------
wald_test(npreg_fit, contrast = riskRatio, null = 1)

## ------------------------------------------------------------------------
DeltaA <- rbinom(n, 1, plogis(2 + W$W1))
DeltaY <- rbinom(n, 1, plogis(2 + W$W2 + A))
A[DeltaA == 0] <- NA
Y[DeltaY == 0] <- NA

## ---- cache = TRUE-------------------------------------------------------
glm_fit_wNAs <- drtmle(W = W, A = A, Y = Y, stratify = FALSE,
                       glm_g = list(DeltaA = "W1 + W2", A = "W1 + W2",
                                    DeltaY = "W1 + W2 + A"),
                       glm_Q = "W1 + W2*A", glm_Qr = "gn",
                       glm_gr = "Qn", family = binomial())
glm_fit_wNAs

## ---- cache = TRUE-------------------------------------------------------
mixed_fit_wNAs <- drtmle(W = W, A = A, Y = Y, stratify = FALSE,
                         SL_g = list(DeltaA = "SL.glm", A = "SL.npreg",
                         DeltaY = c("SL.glm", "SL.mean", "SL.earth")),
                         glm_Q = "W1 + W2*A", glm_Qr = "gn",
                         glm_gr = "Qn", family = binomial())
mixed_fit_wNAs

## ---- echo = FALSE, message = FALSE--------------------------------------
library(drtmle)
set.seed(1234)
N <- 1e6
W <- data.frame(W1 = runif(N), W2 = rbinom(N, 1, 0.5))
Y2 <- rbinom(N, 1, plogis(W$W1 + 2*W$W2))
Y1 <- rbinom(N, 1, plogis(W$W1 + W$W2))
Y0 <- rbinom(N, 1, plogis(W$W1))
EY1 <- mean(Y1)
EY0 <- mean(Y0)
EY2 <- mean(Y2)

## ---- cache = TRUE, message = FALSE--------------------------------------
set.seed(1234)
n <- 200
W <- data.frame(W1 = runif(n), W2 = rbinom(n, 1, 0.5))
A <- rbinom(n, 2, plogis(W$W1 + W$W2))
Y <- rbinom(n, 1, plogis(W$W1 + W$W2*A))

## ---- cache = TRUE, message = FALSE--------------------------------------
glm_fit_multiA <- drtmle(W = W, A = A, Y = Y, stratify = FALSE,
                         glm_g = "W1 + W2", glm_Q = "W1 + W2*A",
                         glm_gr = "Qn", glm_Qr = "gn",
                         family = binomial(), a_0 = c(0,1,2))
glm_fit_multiA

## ---- cache = TRUE-------------------------------------------------------
ci(glm_fit_multiA)
wald_test(glm_fit_multiA, null = c(0.4, 0.5, 0.6))

## ---- cache = TRUE-------------------------------------------------------
ci(glm_fit_multiA, contrast = c(-1, 1, 0))

## ---- cache = TRUE-------------------------------------------------------
ci(glm_fit_multiA, contrast = c(-1, 0, 1))

## ---- cache = TRUE-------------------------------------------------------
riskRatio_1v0 <- list(f = function(eff){ log(eff) },
                      f_inv = function(eff){ exp(eff) },
                      h = function(est){ est[2]/est[1] },
                      fh_grad =  function(est){ c(1/est[2], -1/est[1], 0) })
ci(glm_fit_multiA, contrast = riskRatio_1v0)

## ---- cache = TRUE-------------------------------------------------------
riskRatio_2v0 <- list(f = function(eff){ log(eff) },
                      f_inv = function(eff){ exp(eff) },
                      h = function(est){ est[3]/est[1] },
                      fh_grad =  function(est){ c(0, -1/est[1], 1/est[3]) })
ci(glm_fit_multiA, contrast = riskRatio_2v0)

## ---- cache = TRUE-------------------------------------------------------
cv_sl_fit <- drtmle(W = W, A = A, Y = Y, family = binomial(),
                    SL_g = c("SL.glm", "SL.glm.interaction", "SL.earth"),
                    SL_Q = c("SL.glm", "SL.glm.interaction", "SL.earth"),
                    SL_gr = c("SL.glm", "SL.earth", "SL.mean"),
                    SL_Qr = c("SL.glm", "SL.earth", "SL.mean"),
                    stratify = FALSE, cvFolds = 2, a_0 = c(0, 1, 2))
cv_sl_fit

## ---- cache = TRUE, message = FALSE--------------------------------------
library(future.batchtools)
library(parallel)
cl <- makeCluster(2, type = "SOCK")
plan(cluster, workers = cl)
clusterEvalQ(cl, library("SuperLearner"))
pcv_sl_fit <- drtmle(W = W, A = A, Y = Y, family = binomial(),
                     SL_g = c("SL.glm", "SL.glm.interaction","SL.earth"),
                     SL_Q = c("SL.glm", "SL.glm.interaction","SL.earth"),
                     SL_gr = c("SL.glm", "SL.earth", "SL.mean"),
                     SL_Qr = c("SL.glm", "SL.earth", "SL.mean"),
                     stratify = FALSE, a_0 = c(0,1,2),
                     cvFolds = 2, parallel = TRUE)
pcv_sl_fit

## ---- cache = TRUE-------------------------------------------------------
npreg_iptw_multiA <- adaptive_iptw(W = W, A = A, Y = Y, stratify = FALSE,
                                   SL_g = "SL.npreg", SL_Qr = "SL.npreg",
                                   family = binomial(), a_0 = c(0, 1, 2))
npreg_iptw_multiA

## ------------------------------------------------------------------------
sessionInfo()

