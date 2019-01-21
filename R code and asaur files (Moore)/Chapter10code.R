# # # # # # # # # # # # # # # # # # #
# Applied Survival Analysis Using R
# Dirk F. Moore
# Springer (2016)
# # # # # # # # # # # # # # # # # # #

# Chapter 10

# install.packages("asaur")   # must be done once

# # # # # # # # # # # # # #
# Section 10.3.1   assessing the Weibull distribution
# # # # # # # # # # # # # #

library(asaur)
timeMonths <- gastricXelox$timeWeeks*7/30.25 
delta <- gastricXelox$delta
library(survival) 
result.km <- survfit(Surv(timeMonths, delta) ~ 1)

survEst <- result.km$surv 
survTime <- result.km$time
logLogSurvEst <- log(-log(survEst)) 
logSurvTime <- log(survTime)

# Fig 10.1
plot(logLogSurvEst ~ logSurvTime) 
result.lm <- lm(logLogSurvEst ~ logSurvTime)
abline(result.lm)

attach(pharmacoSmoking) 
ttr[ttr == 0]  <- 0.5 
result.surv <- survfit(Surv(ttr, relapse) ~ 1) 
survEst <- result.surv$surv
survTime <- result.surv$time 

logLogSurvEst <- log(-log(survEst))
logSurvTime <- log(survTime)
result.lm <- lm(logLogSurvEst ~ logSurvTime) 
result.lm

# Fig 10.2
plot(logLogSurvEst ~ logSurvTime) 
abline(result.lm)


# # # # # # # # # # # # # #
# Section 10.3.2   Weibull maximum likelihood
# # # # # # # # # # # # # #

logLikWeib <- function(par, tt, status) {
   mu <- par[1]
   sigma <- par[2]
   lambda.p <- exp(-mu)
   alpha.p <- 1/sigma

   dd <- sum(status)
   sum.t <- sum(status*log(tt))
   sum.t.alpha <- sum(tt^alpha.p)

   term.1 <- dd*log(alpha.p) + alpha.p*dd*log(lambda.p)
   term.2 <- (alpha.p - 1)*sum.t
   term.3 <- (lambda.p^alpha.p)*sum.t.alpha
   result <- term.1 + term.2 - term.3
   result    
   }
   
result <- optim(par=c(4.568, 2.280), fn=logLikWeib, method="L-BFGS-B",
    lower=c(0.001, 0.01), upper=c(5, 5),
    control=list(fnscale = -1),
    tt=ttr, status=relapse)  
result$par

result.survreg.0 <- survreg(Surv(ttr, relapse) ~ 1, dist="weibull")
summary(result.survreg.0)

# # # # # # # # # # # # # #
# Section 10.3.3   Profile Weibull likelihood
# # # # # # # # # # # # # #

logLikWeibProf <- function(par, tt, status) {
   # find log-likelihood for a particular sigma, using mle for mu
   sigma <- par
   alpha.p <- 1/sigma
   dd <- sum(status)
   sum.t <- sum(status*log(tt))
   sum.t.alpha <- sum(tt^alpha.p)
   lambda.p <- (dd/sum.t.alpha)^(1/alpha.p)

   term.1 <- dd*log(alpha.p) + alpha.p*dd*log(lambda.p)
   term.2 <- (alpha.p - 1)*sum.t
   term.3 <- (lambda.p^alpha.p)*sum.t.alpha
   result <- term.1 + term.2 - term.3
   result 
   }
   
resultProf <- optim(par=c(2.280), fn=logLikWeibProf, method="L-BFGS-B",
     lower=c(0.01), upper=c(5), control=list(fnscale = -1),
     tt=ttr, status=relapse)  
sigma.hat <- resultProf$par 
sigma.hat  

dd <- sum(relapse)
sigma <- resultProf$par
alpha.p <- 1/sigma.hat
sum.t.alpha <- sum(ttr^alpha.p)
lambda.p <- (dd/sum.t.alpha)^(1/alpha.p)
mu.hat <- -log(lambda.p)
mu.hat 

sigma.list <- (100:500)/100
n.list <- length(sigma.list)
logLik.list <- rep(NA, n.list)
for (i in 1:n.list) {
  logLik.list[i] <- logLikWeibProf(par=sigma.list[i], ttr, relapse)
  } 

# Fig. 10.3
plot(logLik.list ~ sigma.list, type="l", xlab="sigma", 
   ylab="profile log-likelihood") 
abline(v=sigma.hat, col="gray")

# # # # # # # # # # # # # #
# Section 10.3.4  Selecting a Weibull distribution
# # # # # # # # # # # # # #

result.surv <- survfit(Surv(ttr, relapse) ~ 1, subset={grp =="patchOnly"})
result.summ <- summary(result.surv, time=c(28, 84))
t.vec <- result.summ$time
s.vec <- result.summ$surv
data.frame(t.vec, s.vec)

install.packages("Hmisc")   # must be done once
library(Hmisc)
pharmWeib <- Weibull2(t.vec, s.vec)
t.vals <- 1:200
s.vals <- pharmWeib(t.vals)

model.pharm.weib.basic <- survreg(Surv(ttr, relapse) ~ 1,
     dist="weibull", subset={grp =="patchOnly"} )
mu.hat <- model.pharm.weib.basic$coefficients
sigma.hat <- model.pharm.weib.basic$scale
lambda.hat <- exp(-mu.hat)      # " 1 / scale"
alpha.hat <- 1/sigma.hat        # "shape"
s.mle.vals <- 1 - pweibull(t.vals, shape=alpha.hat, scale=1/lambda.hat)

model.pharm.weib.basic <- survreg(Surv(ttr, relapse) ~ 1,
     dist="weibull", subset={grp =="patchOnly"} )
mu.hat <- model.pharm.weib.basic$coefficients
sigma.hat <- model.pharm.weib.basic$scale
lambda.hat <- exp(-mu.hat)      # " 1 / scale"
alpha.hat <- 1/sigma.hat        # "shape"
s.mle.vals <- 1 - pweibull(t.vals, shape=alpha.hat, scale=1/lambda.hat)

# Fig. 10.4

plot(result.surv, conf.int=F, lwd=2, xlab="Days to relapse",
   ylab="Survival probability", cex.lab=1.3, cex.axis=1.3)
lines(s.mle.vals ~ t.vals, col="blue", lwd=2)
lines(s.vals ~ t.vals, col="red", lwd=2)
points(t.vec, s.vec, col="red", pch=16, cex=1.5)
legend("topright", legend=c("Maximum likelihood", "Match at 24 and 84 days"),
  col=c("blue", "red"), lwd=2)

# # # # # # # # # # # # # #
# Section 10.3.6  accelerated failure times and proportional hazards
# # # # # # # # # # # # # #

result.survreg.grp <- survreg(Surv(ttr, relapse) ~ grp, dist="weibull")
summary(result.survreg.grp)

result.coxph.grp <- coxph(Surv(ttr, relapse) ~ grp)
summary(result.coxph.grp) 

mu0.hat <- result.survreg.grp$coef[1]
sigma.hat <- result.survreg.grp$scale 
alpha.hat <- 1/sigma 
lambda0.hat <- exp(-mu0.hat)
tt.vec <- 0:182 
surv0.vec <- 1 - pweibull(tt.vec, shape=alpha.hat, scale=1/lambda0.hat)

gamma.hat <- result.survreg.grp$coef[2]
surv1.vec <- surv0.vec^(exp(-gamma.hat/sigma.hat))
coxph.surv.est <- survfit(result.coxph.grp, 
    newdata=data.frame(list(grp=c("combination","patchOnly"))))
    
# Fig. 10.5

windows(height=5, width=7)
plot(coxph.surv.est, col=c("red", "black"), xlab="Time in days",
  ylab="Survival probability", cex.axis=1.5, cex.lab=1.5, lwd=2)
lines(surv0.vec ~ tt.vec, col="red", lwd=1)
lines(surv1.vec ~ tt.vec, col="black", lwd=1)
legend("topright", legend=c("Cox model combination", "Cox model patch",
       "Weibull model combination", "Weibull model  patch"),
       col=c("red", "black", "red", "black"), lwd=c(2,2,1,1))

# # # # # # # # # # # # # #
# Section 10.3.7  Weibull model with multiple covariates
# # # # # # # # # # # # # #

modelAll2.coxph <- coxph(Surv(ttr, relapse) ~ grp + age + employment)
summary(modelAll2.coxph) 

model.pharm.weib <- survreg(Surv(ttr, relapse) ~ grp + age + employment, dist="weibull")
summary(model.pharm.weib)

weib.coef.all <- model.pharm.weib$coef 
weib.coef <- weib.coef.all[2:5]
weib.coef.ph <- -weib.coef/model.pharm.weib$scale

model.pharm.coxph <- coxph(Surv(ttr, relapse) ~ grp + age + employment)
coxph.coef <- model.pharm.coxph$coef 
data.frame(weib.coef.ph, coxph.coef)

# # # # # # # # # # # # # #
# Section 10.3.8  Weibull model selection and residual analysis
# # # # # # # # # # # # # #

modelAll.pharm.weib <- survreg(Surv(ttr, relapse) ~ grp + gender + 
                   race + employment + yearsSmoking + levelSmoking +
                   age + priorAttempts + longestNoSmoke, dist="weibull") 
model.step.pharm.weib <- step(modelAll.pharm.weib)

resid.deviance <- residuals(model.pharm.weib, type="deviance")

# the following function is from the appendix:
smoothSEcurve <- function(yy, xx) {
  # use after a call to "plot"
  # fit a lowess curve and 95% confidence interval curve
  xx.list <- min(xx) + ((0:100)/100)*(max(xx) - min(xx)) # make list of x values

  # Then fit loess function through the points (xx, yy) at the listed values
  yy.xx <- predict(loess(yy ~ xx), se=T, newdata=data.frame(xx=xx.list))

  lines(yy.xx$fit ~ xx.list, lwd=2)
  lines(yy.xx$fit - qt(0.975, yy.xx$df)*yy.xx$se.fit ~ xx.list, lty=2)
  lines(yy.xx$fit + qt(0.975, yy.xx$df)*yy.xx$se.fit ~ xx.list, lty=2)
  }

# Fig. 10.6
par(mfrow=c(2,2)) 
plot(resid.deviance ~ age) 
smoothSEcurve(resid.deviance, age) 
title("Deviance residuals\nversus age")

plot(resid.deviance ~ grp)   
title("Deviance residuals\nversus treatment group")

plot(resid.deviance ~ employment) 
title("Deviance residuals\nversus employment") 



resid.dfbeta <- residuals(model.pharm.weib, type="dfbeta")
n.obs <- length(ttr)
index.obs <- 1:n.obs 

# Fig. 10.7

windows()
plot(resid.dfbeta[,3] ~ index.obs, type="h",
   xlab="Observation", ylab="Change in coefficient",
   ylim=c(-0.0065, 0.004)) 
abline(h=0)

identify(resid.dfbeta[,3] ~ index.obs)   # use mouse to click points to identify
#  press "Esc" key to exit

# # # # # # # # # # # # # #
# Section 10.4 Other parametric distributions
# # # # # # # # # # # # # #

model.pharm.lognormal <- survreg(Surv(ttr, relapse) ~ grp + age +
    employment, dist="lognormal")
summary(model.pharm.lognormal)

model.pharm.loglogistic <- survreg(Surv(ttr, relapse) ~ grp + 
    age + employment, dist="loglogistic")
summary(model.pharm.loglogistic)





    
