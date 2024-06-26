
library("RTMB")

set.seed(987)
N <- 50
truth_mean <- c(0.1, 0.8)
truth_sd <- c(0.5, 0.8)
truth <- cbind(
  pnorm(rnorm(N, qnorm(truth_mean[1]), truth_sd[1])),
  pnorm(rnorm(N, qnorm(truth_mean[2]), truth_sd[2]))
)
round(truth[1:10,], 2)
tprob <- cbind(truth[,1], (1-truth[,1])*truth[,2], (1-truth[,1])*(1-truth[,2]))
dat <- matrix(NA_real_, N, 3)
for (i in seq_len(N)) {
  dat[i,] <- rmultinom(1, 100, tprob[i,])
}


datlist <- list(datu = dat,
                nsubj = N)

par_multi <- list(
  mu_p1 = 0,
  mu_p2 = 0,
  sd_p1 = 1,
  sd_p2 = 1,
  p1_u = rep(0, N),
  p2_u = rep(0, N)
)


fun_multi <- function(parms) {
  getAll(parms, datlist, warn=FALSE)
  datu <- OBS(datu)

  nll <- 0
  nll <- nll - sum(dnorm(p1_u, mean=mu_p1, sd=sd_p1, log=TRUE))
  nll <- nll - sum(dnorm(p2_u, mean=mu_p2, sd=sd_p2, log=TRUE))

  ## transform to unit range
  p1 <- pnorm(p1_u)
  p2 <- pnorm(p2_u)

  prob1 <- p1
  prob2 <- (1 - p1) * p2
  prob3 <- (1 - p1) * (1 - p2)

  for (i in 1:nrow(datu)) {
      nll <- nll - dmultinom(datu[i,], prob=c(prob1[i], prob2[i], prob3[i]), log=TRUE)
  }
  nll
}
fun_multi(par_multi)

obj_multi <- MakeADFun(fun_multi, par_multi, random = c("p1_u", "p2_u"))
opt_multi <- nlminb(obj_multi$par, obj_multi$fn, obj_multi$gr)
#zapsmall(opt_multi$par)
sdr_multi <- sdreport(obj_multi)
sdr_multi
summary(sdr_multi, getJointPrecision = TRUE)
pnorm(sdr_multi$par.fixed[1:2])
#     mu_p1     mu_p2
# 0.1076055 0.7732125

round(sdr_multi$par.fixed[3:4], 2)
# sd_p1 sd_p2
#  0.54  0.72

#### check Bayesian method

library("TreeBUGS") ## uses JAGS
model <- "
T1 V1 p1
T1 V2 (1-p1)*p2
T1 V3 (1-p1)*(1-p2)
"
df <- as.data.frame(dat)
m1 <- traitMPT(model, df)
summary(m1)
# Group-level medians of MPT parameters (probability scale):
#          Mean    SD  2.5%   50% 97.5% Time-series SE n.eff  Rhat R_95%
# mean_p1 0.111 0.017 0.081 0.109 0.149          0.001   180 1.022 1.069
# mean_p2 0.769 0.035 0.697 0.771 0.832          0.003   119 1.027 1.088

# Standard deviation of latent-trait values (probit scale) across individuals:
#                  Mean    SD  2.5%   50% 97.5% Time-series SE n.eff  Rhat R_95%
# latent_sigma_p1 0.566 0.066 0.448 0.561 0.708          0.001  3170 1.010 1.035
# latent_sigma_p2 0.747 0.086 0.598 0.740 0.933          0.002  1458 1.001 1.003


## fit model w/o random effects

par_single <- list(
  p1 = 0.5,
  p2 = 0.5
)
dat_single <- list(datu = colSums(dat))

fun_single <- function(parms) {
  getAll(parms, dat_single, warn=FALSE)
  datu <- OBS(datu)
  prob1 <- p1
  prob2 <- (1 - p1) * p2
  prob3 <- (1 - p1) * (1 - p2)

  nll <- 0
  nll <- nll - dmultinom(datu, prob=c(prob1, prob2, prob3), log=TRUE)
  nll
}

obj <- MakeADFun(fun_single, par_single)
opt <- nlminb(obj$par, obj$fn)
sdr <- sdreport(obj)
sdr
# sdreport(.) result
#    Estimate  Std. Error
# p1 0.138600 0.004886402
# p2 0.730671 0.006759472
# Maximum gradient component: 0.0005182904


### also works
fun <- function(parms) {
  getAll(parms, dat_single, warn=FALSE)
  datu <- OBS(datu)
  prob <- c(
    p1,
    (1 - p1) * p2,
    (1 - p1) * (1 - p2)
  )

  nll <- 0
  nll <- nll - dmultinom(datu, prob=prob, log=TRUE)
  nll
}
fun(par_single)
obj2 <- MakeADFun(fun, par_single)
nlminb(obj2$par, obj2$fn)
sdr2 <- sdreport(obj2)
sdr2

### does not work:
# Error in prob < 0 : invalid comparison with complex values
fun1 <- function(parms) {
  getAll(parms, dat_single, warn=FALSE)
  datu <- OBS(datu)
  prob <- rep(NA_real_, 3)
  prob[1] <- p1
  prob[2] <- (1 - p1) * p2
  prob[3] <- (1 - p1) * (1 - p2)
  #prob <- as.numeric(prob)

  nll <- 0
  nll <- nll - dmultinom(datu, prob=prob, log=TRUE)
  nll
}

fun1(par_single)

obj <- MakeADFun(fun1, par_single)
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr

