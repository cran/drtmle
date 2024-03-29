## ----global_options, include=FALSE--------------------------------------------
knitr::opts_chunk$set(warning=FALSE, message=FALSE)

## ---- echo = FALSE, message = FALSE-------------------------------------------
library(drtmle)
set.seed(1234)
N <- 1e6
W <- data.frame(W1 = runif(N), W2 = rbinom(N, 1, 0.5))
Y1 <- rbinom(N, 1, plogis(W$W1 + W$W2))
Y0 <- rbinom(N, 1, plogis(W$W1))
EY1 <- mean(Y1)
EY0 <- mean(Y0)

## -----------------------------------------------------------------------------
set.seed(1234)
n <- 200
W <- data.frame(W1 = runif(n), W2 = rbinom(n, 1, 0.5))
A <- rbinom(n, 1, plogis(W$W1 + W$W2))
Y <- rbinom(n, 1, plogis(W$W1 + W$W2*A))

## -----------------------------------------------------------------------------
glm_fit_uni <- drtmle(W = W, A = A, Y = Y, family = binomial(),
                      glm_g = "W1 + W2", glm_Q = "W1 + W2*A",
                      glm_gr = "Qn", glm_Qr = "gn", stratify = FALSE)
glm_fit_uni

glm_fit_biv <- drtmle(W = W, A = A, Y = Y, family = binomial(),
                      glm_g = "W1 + W2", glm_Q = "W1 + W2*A",
                      glm_gr = "Qn", glm_Qr = "gn", stratify = FALSE,
                      reduction = "bivariate")
glm_fit_biv

## ---- echo = FALSE, message = FALSE-------------------------------------------
library(SuperLearner)

## ---- warning = FALSE, message = FALSE, cache = TRUE--------------------------
set.seed(1234)
sl_fit <- drtmle(W = W, A = A, Y = Y, family = binomial(),
                 SL_g = c("SL.glm","SL.mean"),
                 SL_Q = c("SL.glm","SL.mean"),
                 SL_gr = c("SL.glm", "SL.gam"),
                 SL_Qr = c("SL.glm", "SL.gam"),
                 stratify = FALSE)
sl_fit

## -----------------------------------------------------------------------------
SL.npreg

## ---- cache = TRUE, message = FALSE-------------------------------------------
npreg_fit <- drtmle(W = W, A = A, Y = Y, family = binomial(),
                    SL_g = "SL.npreg", SL_Q = "SL.npreg",
                    SL_gr = "SL.npreg", SL_Qr = "SL.npreg",
                    stratify = TRUE)
npreg_fit

## ---- message = FALSE---------------------------------------------------------
# outcome adaptive propensity score fit using glms
adaptg_fit <- drtmle(W = W, A = A, Y = Y, family = binomial(),
                     glm_Q = ".^2", glm_g = "Q1W + Q0W",
                     adapt_g = TRUE)

# outcome adaptive propensity score fit using super learner
adaptg_sl_fit <- drtmle(W = W, A = A, Y = Y, family = binomial(),
                     SL_Q = c("SL.glm", "SL.mean"), 
                     SL_g = c("SL.glm", "SL.mean"),
                     adapt_g = TRUE)

## ---- cache = TRUE, message = FALSE-------------------------------------------
# fit a GLM for outcome regression outside of drtmle
out_glm_fit <- glm(Y ~ . , data = data.frame(A = A, W), family = binomial())

# generate Qn list
Qn10 <- list(
  # first entry is predicted values setting A = 1
  predict(out_glm_fit, newdata = data.frame(A = 1, W), type = "response"),
  # second entry is predicted values setting A = 0
  predict(out_glm_fit, newdata = data.frame(A = 0, W), type = "response")
)

# pass this list to drtmle to avoid internal estimation of 
# outcome regression (note propensity score and reduced-dimension
# regressions are still estimated internally)
out_fit1 <- drtmle(W = W, A = A, Y = Y, family = binomial(),
                  glm_g = ".", glm_Qr = "gn", glm_gr = "Qn", Qn = Qn10,
                  # crucial to set a_0 to match Qn's construction!
                  a_0 = c(1,0))
out_fit1

# to be completely thorough let's re-order Qn10 and re-run
Qn01 <- list(Qn10[[2]], Qn10[[1]])

out_fit2 <- drtmle(W = W, A = A, Y = Y, family = binomial(),
                  glm_g = ".", glm_Qr = "gn", glm_gr = "Qn", 
                  # use re-ordered Qn list
                  Qn = Qn01,
                  # a_0 has to be reordered to match Qn01!
                  a_0 = c(0,1))
# should be the same as out_fit1, but re-ordered
out_fit2

## ---- cache = TRUE, message = FALSE-------------------------------------------
# fit a GLM for propensity score outside of drtmle
out_glm_fit <- glm(A ~ . , data = data.frame(W), family = binomial())

# get P(A = 1 | W) by calling predict
gn1W <- predict(out_glm_fit, newdata = data.frame(W), type = "response")

# generate gn list
gn10 <- list(
  # first entry is P(A = 1 | W)
  gn1W,
  # second entry is P(A = 0 | W) = 1 - P(A = 1 | W)
  1 - gn1W
)

# pass this list to drtmle to avoid internal estimation of 
# propensity score (note reduced-dimension regressions are 
# still estimated internally)
out_fit3 <- drtmle(W = W, A = A, Y = Y, family = binomial(),
                   glm_Qr = "gn", glm_gr = "Qn", 
                   Qn = Qn10, gn = gn10, 
                   # crucial to set a_0 to match Qn and gn's construction!
                   a_0 = c(1,0))
out_fit3

# to be completely thorough let's re-order gn10 and re-run
gn01 <- list(gn10[[2]], gn10[[1]])

out_fit4 <- drtmle(W = W, A = A, Y = Y, family = binomial(),
                   glm_Qr = "gn", glm_gr = "Qn", 
                   # use re-ordered Qn/gn list
                   Qn = Qn01, gn = gn01, 
                   # a_0 has to be reordered to match Qn01!
                   a_0 = c(0,1))
# should be the same as out_fit3, but re-ordered
out_fit4

## ---- show_targeted_se, cache = TRUE------------------------------------------
glm_fit_uni_nontargeted_se <- 
  drtmle(W = W, A = A, Y = Y, family = binomial(),
         glm_g = "W1 + W2", glm_Q = "W1 + W2*A",
         glm_gr = "Qn", glm_Qr = "gn", stratify = FALSE, 
         targeted_se = FALSE)
# compare to original call        
list(
  glm_fit_uni$drtmle$cov, # based on targeted OR
  glm_fit_uni_nontargeted_se$drtmle$cov # untargeted OR
)

## ---- fit_w_full_cv, warning = FALSE, cache = TRUE----------------------------
glm_fit_uni_cv <- drtmle(W = W, A = A, Y = Y, family = binomial(),
                         glm_g = "W1 + W2", glm_Q = "W1 + W2*A",
                         glm_gr = "Qn", glm_Qr = "gn",
                         se_cv = "full", se_cvFolds = 2, stratify = FALSE,
                         targeted_se = FALSE) # needed for this to work
# compare variance to previously obtained
list(glm_fit_uni$drtmle$cov, glm_fit_uni_cv$drtmle$cov)

## ---- demo_pcv, cache = TRUE, warning = FALSE---------------------------------
set.seed(1234)
sl_fit_pcv <- drtmle(W = W, A = A, Y = Y, family = binomial(),
                     SL_g = c("SL.glm","SL.mean"),
                     SL_Q = c("SL.glm","SL.mean"),
                     SL_gr = c("SL.glm", "SL.gam"),
                     SL_Qr = c("SL.glm", "SL.gam"),
                     stratify = FALSE,
                     se_cv = "partial")
# compare estimates
# pt estimates should be same
# variance should be larger for pcv
list(
  sl_fit, # no cv
  sl_fit_pcv # partial cv
)

## -----------------------------------------------------------------------------
ci(npreg_fit)

## -----------------------------------------------------------------------------
ci(npreg_fit, contrast = c(1,-1))

## -----------------------------------------------------------------------------
riskRatio <- list(f = function(eff){ log(eff) },
                  f_inv = function(eff){ exp(eff) },
                  h = function(est){ est[1]/est[2] },
                  fh_grad =  function(est){ c(1/est[1],-1/est[2]) })
ci(npreg_fit, contrast = riskRatio)

## -----------------------------------------------------------------------------
logitMean <- list(f = function(eff){ qlogis(eff) },
                  f_inv = function(eff){ plogis(eff) },
                  h = function(est){ est[1] },
                  fh_grad = function(est){ c(1/(est[1] - est[1]^2), 0) })
ci(npreg_fit, contrast = logitMean)

## -----------------------------------------------------------------------------
wald_test(npreg_fit, null = c(0.5, 0.6))

## -----------------------------------------------------------------------------
wald_test(npreg_fit, contrast = c(1, -1))
wald_test(npreg_fit, contrast = c(1, -1), null = 0.05)

## -----------------------------------------------------------------------------
wald_test(npreg_fit, contrast = riskRatio, null = 1)

## -----------------------------------------------------------------------------
DeltaA <- rbinom(n, 1, plogis(2 + W$W1))
DeltaY <- rbinom(n, 1, plogis(2 + W$W2 + A))
A[DeltaA == 0] <- NA
Y[DeltaY == 0] <- NA

## ---- cache = TRUE------------------------------------------------------------
glm_fit_wNAs <- drtmle(W = W, A = A, Y = Y, stratify = FALSE,
                       glm_g = list(DeltaA = "W1 + W2", A = "W1 + W2",
                                    DeltaY = "W1 + W2 + A"),
                       glm_Q = "W1 + W2*A", glm_Qr = "gn",
                       glm_gr = "Qn", family = binomial())
glm_fit_wNAs

## ---- cache = TRUE------------------------------------------------------------
mixed_fit_wNAs <- drtmle(W = W, A = A, Y = Y, stratify = FALSE,
                         SL_g = list(DeltaA = "SL.glm", A = "SL.npreg",
                         DeltaY = c("SL.glm", "SL.mean", "SL.gam")),
                         glm_Q = "W1 + W2*A", glm_Qr = "gn",
                         glm_gr = "Qn", family = binomial())
mixed_fit_wNAs

## ---- cache = TRUE------------------------------------------------------------
# first regress indicator of missing A on W
fit_DeltaA <- glm(DeltaA ~ . , data = W, family = binomial())

# get estimated propensity for observing A
ps_DeltaA <- predict(fit_DeltaA, type = "response")

# now regress A on W | DeltaA = 1
fit_A <- glm(A[DeltaA == 1] ~ . , data = W[DeltaA == 1, , drop = FALSE], 
             family = binomial())

# get estimated propensity for receiving A = 1
ps_A1 <- predict(fit_A, newdata = W, type = "response")
# propensity for receiving A = 0
ps_A0 <- 1 - ps_A1

# now regress DeltaY on A + W | DeltaA = 1. Here we are pooling over 
# values of A (i.e., this is what drtmle does if stratify = FALSE). 
# We could also fit two regressions, one of DeltaY ~ W | DeltaA = 1, A = 1
# and one of DeltaY ~ W | DeltaA = 1, A = 0 (this is what drtmle does if
# stratify = TRUE). 
fit_DeltaY <- glm(DeltaY[DeltaA == 1] ~ . , 
                  data = data.frame(A = A, W)[DeltaA == 1, ], 
                  family = binomial())

# get estimated propensity for observing outcome if A = 1
ps_DeltaY_A1 <- predict(fit_DeltaY, newdata = data.frame(A = 1, W), 
                        type = "response")

# get estimated propensity for observing outcome if A = 0
ps_DeltaY_A0 <- predict(fit_DeltaY, newdata = data.frame(A = 0, W), 
                        type = "response")

# now combine all results into a single propensity score 
gn <- list(
  # propensity for DeltaA = 1, A = 1, DeltaY = 1
  ps_DeltaA * ps_A1 * ps_DeltaY_A1,
  # propensity for DeltaA = 1, A = 0, DeltaY = 1
  ps_DeltaA * ps_A0 * ps_DeltaY_A0
)

# pass in this gn to drtmle
out_fit_ps <-  drtmle(W = W, A = A, Y = Y, stratify = FALSE,
                      # make sure a_0/gn are ordered properly!
                      gn = gn, a_0 = c(1, 0),
                      glm_Q = "W1 + W2*A", glm_Qr = "gn",
                      glm_gr = "Qn", family = binomial())
out_fit_ps

## ---- echo = FALSE, message = FALSE-------------------------------------------
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

## ---- cache = TRUE, message = FALSE-------------------------------------------
set.seed(1234)
n <- 200
W <- data.frame(W1 = runif(n), W2 = rbinom(n, 1, 0.5))
A <- rbinom(n, 2, plogis(W$W1 + W$W2))
Y <- rbinom(n, 1, plogis(W$W1 + W$W2*A))

## ---- cache = TRUE, message = FALSE-------------------------------------------
glm_fit_multiA <- drtmle(W = W, A = A, Y = Y, stratify = FALSE,
                         glm_g = "W1 + W2", glm_Q = "W1 + W2*A",
                         glm_gr = "Qn", glm_Qr = "gn",
                         family = binomial(), a_0 = c(0,1,2))
glm_fit_multiA

## ---- cache = TRUE------------------------------------------------------------
ci(glm_fit_multiA)
wald_test(glm_fit_multiA, null = c(0.4, 0.5, 0.6))

## ---- cache = TRUE------------------------------------------------------------
ci(glm_fit_multiA, contrast = c(-1, 1, 0))

## ---- cache = TRUE------------------------------------------------------------
ci(glm_fit_multiA, contrast = c(-1, 0, 1))

## ---- cache = TRUE------------------------------------------------------------
riskRatio_1v0 <- list(f = function(eff){ log(eff) },
                      f_inv = function(eff){ exp(eff) },
                      h = function(est){ est[2]/est[1] },
                      fh_grad =  function(est){ c(1/est[2], -1/est[1], 0) })
ci(glm_fit_multiA, contrast = riskRatio_1v0)

## ---- cache = TRUE------------------------------------------------------------
riskRatio_2v0 <- list(f = function(eff){ log(eff) },
                      f_inv = function(eff){ exp(eff) },
                      h = function(est){ est[3]/est[1] },
                      fh_grad =  function(est){ c(0, -1/est[1], 1/est[3]) })
ci(glm_fit_multiA, contrast = riskRatio_2v0)

## ---- cache = TRUE------------------------------------------------------------
cv_sl_fit <- drtmle(W = W, A = A, Y = Y, family = binomial(),
                    SL_g = c("SL.glm", "SL.glm.interaction", "SL.gam"),
                    SL_Q = c("SL.glm", "SL.glm.interaction", "SL.gam"),
                    SL_gr = c("SL.glm", "SL.mean"),
                    SL_Qr = c("SL.glm", "SL.mean"),
                    stratify = FALSE, cvFolds = 2, a_0 = c(0, 1, 2))
cv_sl_fit

## ---- cache = TRUE, message = FALSE, eval = FALSE-----------------------------
#  ## Commented out to avoid build errors
#  # library(future.batchtools)
#  # library(parallel)
#  # cl <- makeCluster(2, type = "SOCK")
#  # plan(cluster, workers = cl)
#  # clusterEvalQ(cl, library("SuperLearner"))
#  # pcv_sl_fit <- drtmle(W = W, A = A, Y = Y, family = binomial(),
#  #                      SL_g = c("SL.glm", "SL.glm.interaction","SL.gam"),
#  #                      SL_Q = c("SL.glm", "SL.glm.interaction","SL.gam"),
#  #                      SL_gr = c("SL.glm", "SL.gam", "SL.mean"),
#  #                      SL_Qr = c("SL.glm", "SL.gam", "SL.mean"),
#  #                      stratify = FALSE, a_0 = c(0,1,2),
#  #                      cvFolds = 2, use_future = TRUE)

## ---- cache = TRUE------------------------------------------------------------
npreg_iptw_multiA <- adaptive_iptw(W = W, A = A, Y = Y, stratify = FALSE,
                                   SL_g = "SL.npreg", SL_Qr = "SL.npreg",
                                   family = binomial(), a_0 = c(0, 1, 2))
npreg_iptw_multiA

## -----------------------------------------------------------------------------
sessionInfo()

