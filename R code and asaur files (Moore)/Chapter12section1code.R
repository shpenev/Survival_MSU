# # # # # # # # # # # # # # # # # # #
# Applied Survival Analysis Using R
# Dirk F. Moore
# Springer (2016)
# # # # # # # # # # # # # # # # # # #

# # # # # # # # # # # # # # # # #
# Section 12.1 piecewise constant hazards
# # # # # # # # # # # # # # # # #

# install.packages("asaur")   # must be done once
library(survival)
library(asaur)

tt <- c(6, 7, 10, 15, 19, 25)
delta <- c(1, 0, 1, 1, 0, 1)
trt <- c(0, 0, 1, 0, 1, 1)
id <- 1:6 
simple <- data.frame(id, tt, delta, trt) 
simple

tau.s <- c(0, 8, 16, 30)

# # # # # # # # # # # #
# the following code, from the text, has been updated
# simple.split.s <- survSplit(data=simple, cut=tau.s, end="tt",
#      start="t0", event="delta", episode="diagGrp")
# # # # # # # # # # # #
      
simple.split.s <- survSplit(Surv(tt, delta) ~ ., data=simple, cut=tau.s, end="tt",
      start="t0", event="delta", episode="diagGrp")
simple.split.s$expo <- simple.split.s$tt - simple.split.s$t0
ord <- order(simple.split.s$id)
simple.split.ord <- simple.split.s[ord,]
simple.split.ord  

result.simple.poisson <- glm(delta ~ -1 + factor(diagGrp) + trt +
  offset(log(expo)), family=poisson, data=simple.split.ord) 
summary(result.simple.poisson)    

simple.tab <- aggregate(simple.split.ord[c("delta", "expo")],
          by=list(treat=simple.split.ord$trt, 
          diagGrp=simple.split.ord$diagGrp), sum)
simple.tab

result.simple.tab.poisson <- glm(delta ~ -1 + factor(diagGrp) + treat +
     offset(log(expo)), family=poisson, data=simple.tab)

# coxph(Surv(tt, status) ~ grp)  # this code, from the text, is incorrect
coxph(Surv(tt, delta) ~ trt)     # here is the corrected code

alpha0.hat <- as.numeric(result.simple.tab.poisson$coef[1:3])
beta.hat <- result.simple.tab.poisson$coef[4] 
alpha1.hat <- alpha0.hat + beta.hat

install.packages("msm")
library(msm)    
tt.vec <- (0:300)/10 
piece.surv.0 <- ppexp(q=tt.vec, rate=exp(alpha0.hat), t=tau.s[1:3],
   lower.tail=F) 
piece.surv.1 <- ppexp(q=tt.vec, rate=exp(alpha1.hat), t=tau.s[1:3],
   lower.tail=F) 

attach(prostateSurvival)
prost.80plus.poor <- prostateSurvival[{{grade == "poor"}  & {ageGroup == "80+"}},]
prost.80plus.poor$status.all <- as.numeric(prost.80plus.poor$status >= 1) 
prost.80plus.poor$T2 <- as.numeric(prost.80plus.poor$stage == "T2")
prost.80plus.poor$id <- 1:nrow(prost.80plus.poor)
head(prost.80plus.poor)
dim(prost.80plus.poor)
tau.s <- (0:5)*24  
tau.s 

# note new argument:  Surv(survTime, status) ~ .    as required by update in survival package

prost.80plus.poor$survTime[prost.80plus.poor$survTime == 0] <- 0.1  # minimum time must be > 0 for survSplit
prost.split.s <- survSplit(Surv(survTime, status.all) ~ ., data=prost.80plus.poor, cut=tau.s,
   end="survTime", start="t0", event="status.all", episode="survGrp")  
   
prost.split.s$expo <- prost.split.s$survTime - prost.split.s$t0

ord <- order(prost.split.s$id)
prost.split.ord <- prost.split.s[ord,]
prost.tab <- aggregate(prost.split.ord[c("status.all", "expo")],
         by=list(T2=prost.split.ord$T2, 
         survGrp=prost.split.ord$survGrp), sum)
prost.tab

result.prost.tab.poisson <- glm(status.all ~ -1 + 
      factor(survGrp) + T2 + offset(log(expo)), 
      family=poisson, data=prost.tab)
summary(result.prost.tab.poisson) 

alpha0.hat <- as.numeric(result.prost.tab.poisson$coef[1:5])
beta.hat <- result.prost.tab.poisson$coef[6] 
alpha1.hat <- alpha0.hat + beta.hat
library(msm)       # Note: "msm" must have been previously installed
tt.vec <- 0:120
piece.surv.0 <- ppexp(q=tt.vec, rate=exp(alpha0.hat), 
   t=tau.s[1:5], lower.tail=F) 
piece.surv.1 <- ppexp(q=tt.vec, rate=exp(alpha1.hat), 
   t=tau.s[1:5], lower.tail=F)
   
# Fig. 12.4
windows(width=7, height=5)
plot(piece.surv.0 ~ tt.vec, type="n", xlab="Time in months",
    ylab="Survival probability", cex.lab=1.5, cex.axis=1.5)
lines(piece.surv.0 ~ tt.vec, lwd=2)
lines(piece.surv.1 ~ tt.vec, lwd=2, lty=2) 
legend("bottomleft", legend=c("T1", "T2"), lty=c(1,2), lwd=2)

summary(coxph(Surv(survTime, status.all) ~ T2, 
      data=prost.80plus.poor)) 


   