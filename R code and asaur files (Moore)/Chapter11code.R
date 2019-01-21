# # # # # # # # # # # # # # # # # # #
# Applied Survival Analysis Using R
# Dirk F. Moore
# Springer (2016)
# # # # # # # # # # # # # # # # # # #

# Chapter 11

# install.packages("asaur")   # must be done once

# # # # # # # # # # # # # # #
# Sec. 11.1 power and sample size for a single arm study
# # # # # # # # # # # # # # #

expLogMeanDeaths <- function(Delta, alpha, pwr) {   
	z.alpha <- qnorm(alpha, lower.tail=F)   
	z.beta <- qnorm(1-pwr, lower.tail=F)   
	num <- (z.alpha + z.beta)^2   
	denom <- (log(Delta))^2   
	dd <- num/denom   
	dd   
  }
expLikeRatio <- function(d, alpha, pwr) {   
	num <- qchisq(alpha, df=(2*d), lower.tail=F)   
	denom <- qchisq(pwr, df=(2*d), lower.tail=F)   
	Delta <- num/denom   
	Delta 
  }
expLRdeaths <- function(Delta, alpha, pwr) {   
	LRD <- function(x, alpha, pwr) 
             expLikeRatio(x, alpha, pwr) - Delta   
	result <- uniroot(f=LRD, lower=1,upper=1000,  
             alpha=alpha, pwr=pwr)    
	result$root 
  }
  
expLRdeaths(1.5, 0.05, 0.8) 
expLogMeanDeaths(1.5, 0.05, 0.8) 

# # # # # # # # # # # # # # #
# Sec. 11.2 probability of death
# # # # # # # # # # # # # # #
 
prob.death <- function(lambda, accrual, followup) {   
  probDeath <- 1 - (1/(accrual*lambda))*
    (exp(-lambda*followup) - exp(-lambda*(accrual + followup)))   
	probDeath   
	}
prob.death(lambda=0.10, accrual=2, followup=3) 

# # # # # # # # # # # # # # #
# Sec. 11.5 prob of death from non-parametric survival curve
# # # # # # # # # # # # # # #

library(survival)
library(asaur)
attach(gastricXelox)
timeMonths <- timeWeeks*7/30.25

result.km <- survfit(Surv(timeMonths, delta) ~ 1, 
     conf.type="log-log")
timesXe <- result.km$time 
survXe <- result.km$surv

accrual <- 12 
followup <- 6 
times.use <- c(followup, timesXe[{timesXe >= followup} &
   {timesXe <= accrual + followup}]) 
surv.use <- summary(result.km, times=times.use)$surv

times.diff <- diff(c(times.use, accrual + followup)) 
pi.rec <- 1 - (1/accrual)*sum(times.diff*surv.use) 
pi.rec 

surv.simpson <- summary(result.km, 
      times=c(followup, accrual/2 + followup, accrual+followup))$surv 
surv.simpson 

pi.simpson <- 1 - (1/6)*(surv.simpson[1] + 4*surv.simpson[2] + surv.simpson[3]) 

# # # # # # # # # # # # # # #
# Sec. 11.6 required number of patients for advanced gastric cancer
# # # # # # # # # # # # # # #

TwoArmDeaths <- function(Delta, p=0.5, alpha=0.025, pwr=0.8) {    	
  z.alpha <- qnorm(alpha, lower.tail=F)    	
  z.beta <- qnorm(1-pwr, lower.tail=F)    	
  num <- (z.alpha + z.beta)^2    	
  denom <- p*(1-p)*(log(Delta))^2    	
  dd <- num/denom    	
  dd   
  }
  
TwoArmDeaths(Delta=2, p=0.5, alpha=0.025, pwr=0.8) 
Delta <- 2 
surv.alt <- surv.use^(1/Delta) 
surv.avg <- 0.5*surv.use + 0.5*surv.alt
pi.exact <- 1 - (1/accrual)*sum(times.diff*surv.avg) 
pi.exact

pi0 <- prob.death(lambda=0.0673, accrual=12, followup=6)  
pi1 <- prob.death(lambda=0.0336, accrual=12, followup=6)  
pi.harmonicMean <- 1/(0.5/pi0 + 0.5/pi1) 
pi.harmonicMean 

pi.avg <- (pi0 + pi1)/2
pi.avg 

# # # # # # # # # # # # # # #
# Sec. 11.7 required number of patients for metastatic colorectal cancer
# # # # # # # # # # # # # # #

TwoArmDeaths(Delta=1.834, p=0.5, alpha=0.025, pwr=0.85) 
pDeathSimpson <- function(aa, ff, S) { 
 #Use Simpson's rule to approximate the probability of death 
 #   assuming uniform accrual  
 probDeath <- 1 - (1/6)*(S[1] + 4*S[2] + S[3])   
 probDeath   
 }
 
aa=2 
ff=2 
So <- c(0.76, 0.59, 0.49) 
psi = 1.834 
Sa <- So^(1/psi) 
Sboth <- 0.5*(So + Sa)

pDeathControl <- pDeathSimpson(aa=2, ff=2, S=So)
pDeathControl 
pDeathTreatment <- pDeathSimpson(aa=2, ff=2, S=Sa)
pDeathTreatment 
pDeathAll <- 0.5*(pDeathControl + pDeathTreatment)
pDeathAll 

# # # # # # # # # # # # # # #
# Sec. 11.8 power simulations
# # # # # # # # # # # # # # #

library(survival)
library(asaur)

attach(prostateSurvival)
prost.66to74.poor <- prostateSurvival[{{grade == "poor"}  & 
    {{ageGroup == "70-74"} | {ageGroup == "66-69"}} & {stage == "T2"}},]
library(survival)
status.prost <- as.numeric(prost.66to74.poor$status == 1)
result.prost.survfit <- survfit(Surv(survTime, status.prost) ~ 1, 
    data=prost.66to74.poor)
summary(result.prost.survfit, time=c(48, 96)) 

install.packages("Hmisc")  # must be done once
library(Hmisc) 
Weib.p <- Weibull2(c(4,8),c(0.931,0.717))

Weib.p 
function (times = NULL, alpha = 0.0033021237632906, gamma = 2.21819823268731) {
    exp(-alpha * (times^gamma)) 
	} 
ff <- Quantile2(Weib.p,hratio=function(x) 0.75) 

# Fig. 11.3
plot(ff, xlim=c(0,8))

rcontrol <- function(n) ff(n, what='control') 
rintervention  <- function(n) ff(n, what='intervention')
rcens <- function(n) runif(n, 5, 8)
spower(rcontrol, rintervention, rcens, nc=1500, ni=1500,  
        test=logrank, nsim=1000, alpha=0.025) 
 


