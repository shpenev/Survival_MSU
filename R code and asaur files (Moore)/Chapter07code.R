# # # # # # # # # # # # # # # # # # #
# Applied Survival Analysis Using R
# Dirk F. Moore
# Springer (2016)
# # # # # # # # # # # # # # # # # # #

# Chapter 7

# install.packages("asaur") # must be run once

# # # # # # # # # # # # # #
# Section 7.1.1 martingale and deviance residuals
# # # # # # # # # # # # # #
library(asaur)
attach(pharmacoSmoking)
priorAttemptsT <- priorAttempts
priorAttemptsT[priorAttempts > 20] <- 20 

library(survival)
result.0.coxph <- coxph(Surv(ttr, relapse) ~ 1)
rr.0 <- residuals(result.0.coxph, type="martingale") 

# Figure 7.1

# First, we need the "smoothSEcurve"  function (also in the Appendix)
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

# Now the plot
par(mfrow=c(3,2))
plot(rr.0 ~ age)
smoothSEcurve(rr.0, age)
title("Martingale residuals\nversus age")

logAge <- log(age)
plot(rr.0 ~ logAge)
smoothSEcurve(rr.0, logAge)
title("Martingale residuals\nversus log age")

plot(rr.0 ~ priorAttemptsT)
smoothSEcurve(rr.0, priorAttemptsT)
title("Martingale residuals versus\nprior attempts")

logPriorAttemptsT <- log(priorAttemptsT + 1)
plot(rr.0 ~ logPriorAttemptsT)
smoothSEcurve(rr.0, logPriorAttemptsT)
title("Martingale residuals versus\nlog prior attempts")

plot(rr.0 ~ longestNoSmoke)
smoothSEcurve(rr.0, longestNoSmoke)
# Note that "\n"  is a "newline" indicator
title("Martingale residuals versus\nlongest period without smoking")

logLongestNoSmoke <- log(longestNoSmoke+1)
plot(rr.0 ~ logLongestNoSmoke)
smoothSEcurve(rr.0, logLongestNoSmoke)
title("Martingale residuals versus\nlog of longest period without smoking") 

# # # # #

result.grp.coxph <- coxph(Surv(ttr, relapse) ~ grp) 
result.step <- step(result.grp.coxph, scope=list(upper=~ grp +
                  gender + race + employment + yearsSmoking + 
                  levelSmoking + age + priorAttemptsT + 
                  logLongestNoSmoke, lower=~grp) )
result.step

# Figure 7.2

rr.final <- residuals(result.step, type="martingale") 
par(mfrow=c(2,2))
plot(rr.final ~ age)
smoothSEcurve(rr.final, age)
title("Martingale residuals\nversus age")

plot(rr.final ~ grp)
title("Martingale residuals\nversus treatment group")

plot(rr.final ~ employment)
title("Martingale residuals\nversus employment") 

# # # # # # # # # # # # # #
# Section 7.1.2 case deletion residuals
# # # # # # # # # # # # # #

result.coxph <- coxph(Surv(ttr, relapse) ~ grp + employment + age) 
coef.all <- result.coxph$coef[4]  
coef.all

n.obs <- length(ttr) 
jkbeta.vec <- rep(NA, n.obs) 
for (i in 1:n.obs) {
  tt.i <- ttr[-i]
  delta.i <- relapse[-i]
  grp.i <- grp[-i]
  employment.i <- employment[-i]
  age.i <- age[-i]
  result.coxph.i <- coxph(Surv(tt.i, delta.i) ~ grp.i + 
        employment.i + age.i)
  coef.i <- result.coxph.i$coef[4]
  jkbeta.vec[i] <- (coef.all - coef.i)
  }

# Fig. 7.3  
windows()
index.obs <- 1:n.obs 
plot(jkbeta.vec ~ index.obs, type="h",
    xlab="Observation", ylab="Change in coefficient for age",
   cex.axis=1.3, cex.lab=1.3) 
abline(h=0)

identify(jkbeta.vec ~ index.obs)  # click on points to get their identifying label.
#   Press the "Esc" button to quit the "identify" mode.

# # # # # # # # # # # # # #
# Section 7.2.1 checking the proportional hazards assumption; log cumulative hazard plots
# # # # # # # # # # # # # #

attach(pancreatic2)
pfs.month <- pfs/30.25
result.surv.LA <- survfit(Surv(pfs.month) ~ stage, subset={stage == "LA"})
time.LA <- result.surv.LA$time
surv.LA <- result.surv.LA$surv
cloglog.LA <- log(-log(surv.LA))
logtime.LA <- log(time.LA)

result.surv.M <- survfit(Surv(pfs.month) ~ stage, subset={stage == "M"})
time.M <- result.surv.M$time
surv.M <- result.surv.M$surv
cloglog.M <- log(-log(surv.M))
logtime.M <- log(time.M)

# Fig 7.4
windows()
plot(cloglog.LA ~ logtime.LA, type="s", col="blue", 
   lwd=2, xlab="Log time", ylab="Complementary log-log survival")
lines(cloglog.M ~ logtime.M, col="red", lwd=2, type="s")
legend("bottomright", legend=c("Locally advanced", "Metastatic"), col=c("blue","red"), lwd=2)

# # # # # # # # # # # # # #
# Section 7.2.2 checking the proportional hazards assumption; Schoenfeld residuals
# # # # # # # # # # # # # #

tt <- c(6, 7, 10, 15, 19, 25) 
delta <- c(1, 0, 1, 1, 0, 1) 
trt <- c(0, 0, 1, 0, 1, 1) 
result.coxph <- coxph(Surv(tt, delta) ~ trt) 
result.coxph$coef 
residuals(result.coxph, type="schoenfeld")

resid.unscaled <- residuals(result.coxph, type="schoenfeld")
resid.scaled <- resid.unscaled*result.coxph$var*sum(delta)

resid.unscaled  
resid.scaled
resid.scaled + result.coxph$coef 

resid.sch <- cox.zph(result.coxph)
resid.sch$y 

result.coxph <- coxph(Surv(pfs.month) ~ stage) 
result.sch.resid <- cox.zph(result.coxph, transform="km")

# Fig 7.5
plot(result.sch.resid) 

result.sch.resid    
cox.zph(result.coxph, transform="rank")
cox.zph(result.coxph, transform="identity")

  