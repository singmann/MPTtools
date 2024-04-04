
library("RTMB")
data( sim, package = "mptmem" )
head( sim )
nrow(sim)

model_text <- "
T1 C11 c*r
T1 C14 c*(1-r)
T1 C12 (1-c)*u*u
T1 C13 (1-c)*u*(1-u)
T1 C13 (1-c)*(1-u)*u
T1 C14 (1-c)*(1-u)*(1-u)"

tmp <- read_mpt(text = model_text, type = "eqn")

mp <- parse_model_df(tmp)
cat(paste0(find.MPT.params(mp), "_u = rep(0, nsubj)", collapse = ",\n"))
cat(paste0(find.MPT.params(mp), "_mu = 0", collapse = ",\n"))
cat(paste0(find.MPT.params(mp), "_sd = 1", collapse = ",\n"))

dat_use <- sim[1:10, attr(mp, "cat_map")[[1]]]
dat_use <- sim[, attr(mp, "cat_map")[[1]]]
nsubj <- nrow(dat_use)

parameters <- list(
  c_sd = 1,
  r_sd = 1,
  u_sd = 1,
  c_mu = 0,
  r_mu = 0,
  u_mu = 0,
  c_u = rep(0, nsubj),
  r_u = rep(0, nsubj),
  u_u = rep(0, nsubj)
)

tmpars <- find.MPT.params(mp)
cat(paste0("nll <- nll - sum(dnorm(",
           tmpars, "_u, mean=",
           tmpars, "_mu, sd=",
           tmpars, "_sd, log=TRUE))",
           collapse = "\n"))
cat(paste0(tmpars, " <- pnorm(", tmpars, "_u)",
           collapse = "\n"))

cat(make_model_equations(tmp))

dlist <- list(
  nsubj = nrow(dat_use),
  dat = as.matrix(dat_use)
)
str(dlist)
#dput(dlist)

fmulti <- function(parms) {
  getAll(dlist, parms, warn=FALSE)
  dat <- OBS(dat)

  nll <- 0
  nll <- nll - sum(dnorm(c_u, mean=c_mu, sd=c_sd, log=TRUE))
  nll <- nll - sum(dnorm(r_u, mean=r_mu, sd=r_sd, log=TRUE))
  nll <- nll - sum(dnorm(u_u, mean=u_mu, sd=u_sd, log=TRUE))

  c <- pnorm(c_u)
  r <- pnorm(r_u)
  u <- pnorm(u_u)

  prob1 <- c*r
  prob2 <- (1-c)*u*u
  prob3 <- (1-c)*u*(1-u) + (1-c)*(1-u)*u
  prob4 <- c*(1-r) + (1-c)*(1-u)*(1-u)
  # subj <- 1
  # nll <- nll - dmultinom(dat, prob = c(prob1[subj],
  #                                          prob2[subj],
  #                                          prob3[subj],
  #                                          prob4[subj]), log = TRUE)

  #prob1 <- matrix(0, nsubj, 4)
  for (subj in seq_len(nsubj)) {
    nll <- nll - dmultinom(dat[subj,], prob = c(prob1[subj],
                                                    prob2[subj],
                                                    prob3[subj],
                                                    prob4[subj]), log = TRUE)
  }
  nll
}
fmulti(parameters)
obj <- MakeADFun(fmulti, parameters, random = paste0(tmpars, "_u"))
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr
# $par
#      c_sd      r_sd      u_sd      c_mu      r_mu      u_mu
# 0.3949489 0.3728911 0.3168043 0.3392110 0.3533142 0.1464151
#
# $objective
# [1] 1040.912

pnorm(sdr$par.fixed[4:6])
#      c_mu      r_mu      u_mu
# 0.6012074 0.6887813 0.5035756

library( mptmem )

df <- sim

m <- '

%Within

#Pairs:
T1:C11~c*r
T1:C14~c*(1-r)
T1:C12~(1-c)*u*u
T1:C13~(1-c)*u*(1-u)
T1:C13~(1-c)*(1-u)*u
T1:C14~(1-c)*(1-u)*(1-u)
#Singletons:
T2:C21~a
T2:C22~(1-a)

%Person
c~~c
r~~r
u~~u
a~~a
'

fit <- mptmem( m, data = df, method = "AGH", type_mpt = "logit" )
summary( fit )


fmulti <- function(parms) {
  getAll(dlist, parms, warn=FALSE)
  dat <- OBS(dat)

  nll <- 0
  nll <- nll - sum(dnorm(c_u, mean=c_mu, sd=c_sd, log=TRUE))
  nll <- nll - sum(dnorm(r_u, mean=r_mu, sd=r_sd, log=TRUE))
  nll <- nll - sum(dnorm(u_u, mean=u_mu, sd=u_sd, log=TRUE))

  c <- pnorm(c_u)
  r <- pnorm(r_u)
  u <- pnorm(u_u)

  #prob1 <- matrix(0, nsubj, 4)
  for (subj in seq_len(nsubj)) {
    prob1 <- c[subj]*r[subj]
    prob2 <- (1-c[subj])*u[subj]*u[subj]
    prob3 <- (1-c[subj])*u[subj]*(1-u[subj]) + (1-c[subj])*(1-u[subj])*u[subj]
    prob4 <- c[subj]*(1-r[subj]) + (1-c[subj])*(1-u[subj])*(1-u[subj])
    nll <- nll - dmultinom(x = dat[subj,], prob = c(prob1, prob2, prob3, prob4), log = TRUE)
  }
  nll
}

fmulti <- function(parms) {
  getAll(dlist, parms, warn=FALSE)

  nll <- 0
  nll <- nll - sum(dnorm(c_u, mean=c_mu, sd=c_sd, log=TRUE))
  nll <- nll - sum(dnorm(r_u, mean=r_mu, sd=r_sd, log=TRUE))
  nll <- nll - sum(dnorm(u_u, mean=u_mu, sd=u_sd, log=TRUE))

  c <- pnorm(c_u)
  r <- pnorm(r_u)
  u <- pnorm(u_u)

  prob1 <- matrix(0, nsubj, 4)
  for (subj in seq_len(nsubj)) {
    prob1[subj, 1] <- c[subj]*r[subj]
    prob1[subj, 2] <- (1-c[subj])*u[subj]*u[subj]
    prob1[subj, 3] <- (1-c[subj])*u[subj]*(1-u[subj]) + (1-c[subj])*(1-u[subj])*u[subj]
    prob1[subj, 4] <- c[subj]*(1-r[subj]) + (1-c[subj])*(1-u[subj])*(1-u[subj])
    nll <- nll - dmultinom(x = dat[subj,], prob = prob1[subj,], log = TRUE)
  }
  nll
}

fmulti(parameters)
obj <- MakeADFun(fmulti, parameters, random = paste0(tmpars, "_u"))
