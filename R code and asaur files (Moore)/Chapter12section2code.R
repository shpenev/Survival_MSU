# # # # # # # # # # # # # # # # # # #
# Applied Survival Analysis Using R
# Dirk F. Moore
# Springer (2016)
# # # # # # # # # # # # # # # # # # #

# # # # # # # # # # # # # # # # #
# Chapter 12 Sec 12.2  Interval censoring
# # # # # # # # # # # # # # # # #

install.packages("ssym")  # must be done once
library(ssym)

# Baboon descent times example
data(Baboons)
Baboons[c(1,39,71,101, 150),]

Baboons <- within(Baboons, {   
  delta <- rep(0, length(cs))   
  delta[cs == 0] <- 1
  tt.L <- t   
  tt.R <- t
  tt.L[cs == 1] <- 0.1 })
Baboons[c(1,39,71,101, 150),]

# install.packages("Icens")      # this package is now in Bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite("Icens")

install.packages("interval")   # must be done once

library(Icens)
library(interval)       
result.icfit <- icfit(Surv(time=tt.L, time2=tt.R, 
     type="interval2") ~ 1, conf.int=T, data=Baboons) 

# Fig. 12.5  non-parametric survival curve only
plot(result.icfit, XLAB="Time in hours", 
    YLAB="Survival probability", estpar=list(col="blue", lwd=2), 
    cipar=list(col="blue", lty=2))

baboon.survreg <- survreg(Surv(time=tt.L, time2=tt.R, 
     type="interval2") ~ 1, dist="weibull", data=Baboons)

ones <- rep(1,nrow(Baboons))
baboon.survreg <- survreg(
     Surv(time=tt.L, time2=tt.R, type="interval2") ~ ones,
     dist="weibull", data=Baboons)
pct <- 1:999/1000   
ptime <- predict(baboon.survreg, type='quantile', 
     newdata=data.frame(ones=1), p=pct, se=TRUE)  

# Fig. 12.5   # add Weibull estimate of survival function
matlines(cbind(ptime$fit, ptime$fit + 2*ptime$se.fit,
        ptime$fit- 2*ptime$se.fit), 1-pct,
        xlab="Hours", ylab="Survival", type='l', lty=c(1,2,2),
        lwd=c(2,1,1), xlim=c(0,20), col="red")

# reverse survival estimate and plot (plot not shown in text)
result.surv.reverse <- survfit(Surv(-t, delta) ~ 1, conf.int=T, 
         data=Baboons, conf.type="log-log") 
plot(result.surv.reverse, xlim=c(0, -18), fun="event")

# Breast Cosmesis Study

# install.packages("interval")   # must be done once
library(interval)
data(bcos)
bcos[c(1,33, 47, 62, 90),]

icout <- icfit(Surv(left,right,type="interval2")~treatment, 
           data=bcos, conf.int=F)  
# Fig. 12.6
plot(icout, XLAB="Time in months", YLAB="Survival probability",
        COL=c("lightblue", "pink"), LEGEND=F,
        estpar=list(col=c("blue", "red"), lwd=2, lty=1))
legend("bottomleft", 
        legend=c("Radiation alone", "Radiation and chemo"),
        col=c("blue","red"), lwd=2)
        
# Now fit a Weibull distribution
bcos <- within(bcos, { 
          left.alt <- left 
          left.alt[left == 0] <- 0.1 
          right.alt <- right 
          right.alt[is.infinite(right)] <- 65})
bcos[c(1,33, 47, 62, 90),]

# add Weibull fitted curve to Fig. 12.6

bcos.survreg <- 
  survreg(Surv(left.alt, right.alt, type="interval2") ~ treatment,
  dist="weibull", data=bcos)
pct <- 1:999/1000   
ptime <- predict(bcos.survreg, type='quantile',
       newdata=data.frame(treatment=c("Rad", "RadChem")),
       p=pct, se=F)                       
lines(ptime[1,], 1-pct, xlab="Hours", ylab="Survival", type='l',
          lty=c(1,2,2), lwd=c(2,1,1), xlim=c(0,20), col="blue") 
lines(ptime[2,], 1-pct, xlab="Hours", ylab="Survival", type='l',
          lty=c(1,2,2), lwd=c(2,1,1), xlim=c(0,20), col="red")

        
