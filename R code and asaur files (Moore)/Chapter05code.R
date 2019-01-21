# # # # # # # # # # # # # # # # # # #
# Applied Survival Analysis Using R
# Dirk F. Moore
# Springer (2016)
# # # # # # # # # # # # # # # # # # #

# Chapter 5

# install.packages("asaur")  # this must be done once
library(survival)

# # # # # # # # # #
# Section 5.2, simple example, 
#  data from Section 4.1
# # # # # # # # # #

tt <- c(6, 7, 10, 15, 19, 25)
delta <- c(1, 0, 1, 1, 0, 1)
trt <- c(0, 0, 1, 0, 1, 1)

plsimple <- function(beta) {
  psi <- exp(beta)
  result <- log(psi) - log(3*psi + 3) - log(3*psi + 1) - log(2*psi + 1)
  result
  }

beta.vec <- ((-200):50)/50
plikel.vec <- plsimple(beta.vec)

# # # # # # # # # #
# Figure 5.2.1
# # # # # # # # # #

plot(plikel.vec ~ beta.vec, type="l", xlab="beta", ylab="log partial likelihood",
   cex.axis=1.5, cex.lab=1.5, lwd=2, col="black", xlim=c(-4, 1.4), ylim=c(-5.5, -3.5))

result <- optim(par=0, fn = plsimple, method = "L-BFGS-B",
                control=list(fnscale = -1), lower = -3, upper = 1)
result$par
abline(v=result$par, lty=2)
abline(v=0, lty=2)

# # # # # # # # #
# Section 5.3.3, the likelihood ratio test
# # # # # # # # #

betahat <- result$par
plmax <- plsimple(betahat)
plmax

library(survival)
status <- delta
grp <- trt
result.cox <- coxph(Surv(tt, status) ~ grp)
summary(result.cox)

install.packages("numDeriv")    # must be done once
library(numDeriv)   # must be downloaded from CRAN and installed
grad(func=plsimple, x=0)
hessian(func=plsimple, x=result$par)
1/sqrt(-hessian(func=plsimple, x=result$par))  # standard error

pchisq(1.274, df=1, lower.tail=F) # score test

wald <- result$par*sqrt(-hessian(func=plsimple, x=result$par))
wald                  # Wald test
result.cox$wald.test  # compare to model output

betahat <- result$par
2*(plsimple(betahat) - plsimple(0))   # likelihood ratio test

# # # # # # # # # # # # #
# Section 5.6 tied survival times
# # # # # # # # # # # # #

tt <- c(7, 6, 6, 5, 2, 4, 4, 1, 3, 1)
status <- c(0, 1, 0, 0, 1, 1, 1, 1, 0, 1)
grp <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1)

library(survival)
result.coxph <- coxph(Surv(tt, status) ~ grp, ties="exact")
summary(result.coxph)

loglikContinuous <- function(b) {
  result <- 3*b + log(exp(b) + 9) - log(4*exp(b) + 6) - log(3*exp(b) + 6) -
    log(2*exp(b) + 6) - log(exp(b) + 5) - log(exp(b) + 4)
  result
  }
loglikDiscrete <- function(b) {
 
  resultA <- exp(2*b)/(6*exp(2*b) + 24*exp(b) + 15)
  resultB <- 1/(6 + 2*exp(b))
  resultC <- exp(b)/(10+5*exp(b))

  result <- log(resultA) + log(resultB) + log(resultC)
  result
  }

result.optim.continuous <- optim(par=1.4, fn=loglikContinuous, method="BFGS",
  control=list(fnscale = -1) )
result.optim.discrete <- optim(par=1.4, fn=loglikDiscrete, method="BFGS",
  control=list(fnscale = -1) )
result.optim.continuous$par
# [1] 1.838603
result.optim.discrete$par
# [1] 1.856768
result.coxph$coef     # same as "discrete" method

# # # # # # # # # # # # # #
# Dection 5.7 Left truncation
# # # # # # # # # # # # # #


tt <- c(6, 7, 10, 15, 19, 25)
status <- c(1, 0, 1, 1, 0, 1)
grp <- c(0, 0, 1, 0, 1, 1)
backTime <- c(-3, -11, -3, -7, -10, -5)

library(survival)
coxph(Surv(tt, status) ~ grp)


tm.enter <- -backTime
tm.exit <- tt - backTime
data.frame(tm.enter, tm.exit, status, grp)

coxph(Surv(tm.enter, tm.exit, status) ~ grp)

# # # # # # # # # # # # #
# Channing House data
# # # # # # # # # # # # #

library(asaur)

# add entryYears and exitYears to data set, as in Section 3.5
ChanningHouse <- within(ChanningHouse, {
   entryYears <- entry/12
   exitYears <- exit/12})

channing68 <- ChanningHouse[ChanningHouse$exitYears >= 68,]

coxph(Surv(entryYears, exitYears, cens) ~ sex, data=channing68)

  





