
library("RTMB")

set.seed(987)
truth <- c(0.1, 0.8)
tprob <- c(truth[1], (1-truth[1])*truth[2], (1-truth[1])*(1-truth[2]))
dat <- rmultinom(1, 100, tprob)


par_single <- list(
  p1 = 0.5,
  p2 = 0.5
)
dat_single <- list(datu = dat[,1])

## works
fun1 <- function(parms) {
  getAll(parms, dat_single, warn=FALSE)
  datu <- OBS(datu)
  prob1 <- p1
  prob2 <- (1 - p1) * p2
  prob3 <- (1 - p1) * (1 - p2)

  nll <- 0
  nll <- nll - dmultinom(datu, prob=c(prob1, prob2, prob3), log=TRUE)
  nll
}

obj <- MakeADFun(fun1, par_single)
opt <- nlminb(obj$par, obj$fn)
sdr <- sdreport(obj)
sdr
# sdreport(.) result
#     Estimate Std. Error
# p1 0.1000002 0.02999870
# p2 0.6999998 0.04830441
# Maximum gradient component: 0.0002285503


### works now:
fun2 <- function(parms) {
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
obj <- MakeADFun(fun2, par_single)
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr

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


