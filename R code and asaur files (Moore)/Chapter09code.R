# # # # # # # # # # # # # # # # # # #
# Applied Survival Analysis Using R
# Dirk F. Moore
# Springer (2016)
# # # # # # # # # # # # # # # # # # #

# Chapter 9

# install.packages("asaur")   # must be run once
install.packages("timereg")   # must be run once, for the "diabetes" data
library(asaur)

# Section9.1      clustered survival and frailty models
ashkenazi[ashkenazi$famID %in% c(1, 9, 94), ]

library(timereg)   # for access to the "diabetes"  data
data(diabetes)
head(diabetes)

# Section 9.1.3     family-based clusters in "ashkenazi" data

result.coxph <- coxph(Surv(age, brcancer) ~ mutant, data=ashkenazi)
summary(result.coxph)
result.coxph$loglik 

result.coxph.cluster <-  coxph(Surv(age, brcancer) ~ mutant +  cluster(famID), data=ashkenazi)
summary(result.coxph.cluster) 

result.coxph.frail <-  coxph(Surv(age, brcancer) ~ mutant + frailty(famID), data=ashkenazi)
summary(result.coxph.frail)  

install.packages("coxme")
library(coxme)
result.coxme <- coxme(Surv(age, brcancer) ~ mutant + (1|famID), data=ashkenazi)
summary(result.coxme)

pchisq(4.246,1,lower.tail=F) 

# Section 9.1.4 within-person pairing of eye observations
result.coxme <- coxme(Surv(time, status) ~ treat + 
   as.factor(adult) + treat*as.factor(adult) + (1 | id), data=diabetes)
summary(result.coxme) 

# # # # # # # # # # # # # # # # # # #
# 9.2.1 competing risks   # prostateSurvival is in "asaur" package
# # # # # # # # # # # # # # # # # # #

library(survival)
library(asaur)
prostateSurvival <- within(prostateSurvival, {
   status.prost <- as.numeric({status == 1}) 
   status.other <- as.numeric({status == 2})}) 
attach(prostateSurvival)
prostateSurvival.highrisk <- prostateSurvival[{{grade == "poor"} & 
   {stage=="T2"} & {ageGroup == "80+"}},]
head(prostateSurvival.highrisk) 

result.prostate.km <- survfit(Surv(survTime, event=status.prost) ~ 1,
    data=prostateSurvival.highrisk)
result.other.km <- survfit(Surv(survTime, event=status.other) ~ 1,
    data=prostateSurvival.highrisk) 

surv.other.km <- result.other.km$surv 
time.km <- result.other.km$time/12 
surv.prost.km <- result.prostate.km$surv
cumDist.prost.km <- 1 - surv.prost.km

# Figure 9.1   (simple)
plot(cumDist.prost.km ~ time.km, type="s", ylim=c(0,1), lwd=2, 
    xlab="Years from prostate cancer diagnosis",  col="blue") 
lines(surv.other.km ~ time.km, type="s", col="green", lwd=2) 

# Figure 9.1 (more complete)
windows(height=5, width=7)
layout(matrix(c(2,1,3), ncol=3),
  heights=c(5,5,5), widths=c(0.5, 5, 0.5))
layout.show(3)
par(mar=c(5,2,2,2), cex.axis=1.65, cex.lab=1.65)

#par(mar=c(5, 4, 4, 4) + 0.1)
plot(cumDist.prost.km ~ time.km, type="s", ylim=c(0,1), lwd=2, xlab="Years from prostate cancer diagnosis",
   axes=F, col="blue")
lines(surv.other.km ~ time.km, type="s", col="green", lwd=2)
axis(1)
axis(2, at=seq(0,1,0.2), las=1)
#axis(3)
axis(4, at=seq(0,1,0.2), labels=seq(1,0,-0.2), las=1)
box()
text(20/12, 0.25, label="Death from\nprostate cancer", cex=1.65)
text(20/12, 0.73, label="Death from\nother causes", cex=1.65)

x <- c(0,5)
y <- c(0,0.25)
par(mar=c(0,0,0,0))
plot(x ~ y,type="n",axes=F)
text(x=0.1,y=2.5,"Probability of death from prostate cancer",srt=90,cex=1.65)
plot(x ~ y,type="n",axes=F)
text(x=0.15,y=2.5,"Probability of death from other causes",srt=-90,cex=1.65)


# # # # # # # # # # # #
# Section 9.2.2 cause-specific hazards and cumulative incidence functions
# # # # # # # # # # # #

tt <- c(2,7,5,3,4,6)
status <- c(1,2,1,2,0,0) 
status.any <- as.numeric(status >= 1) 
result.any <- survfit(Surv(tt, status.any) ~ 1) 
result.any$surv 

install.packages("mstate") # must be done once
library(mstate)

# # # # # # # # # # # # # # # # # # # # # # # # # # # #  ! !
# # # # # # # # # # # # # # # # # # # # # # # # # # # #  ! !
# Note: The following plots render correctly with R version 3.2
# With R version 3.3, released in May, 2016, the cumulative incidence
#  estimates may be given for only a few time points, so that the
#  cumulative incidence plots are too coarse
#
# Thus, it is best to use 
#  survfit(Surv(survTime, status, type="mstate") ~ 1, data=prostateSurvival.highrisk)
# as described below
# # # # # # # # # # # # # # # # # # # # # # # # # # # #  ! !
# # # # # # # # # # # # # # # # # # # # # # # # # # # #  ! !


#ci <- Cuminc(time=tt, status=status)
#ci

survfit.mstate.simple  <- survfit(Surv(tt, status, type="mstate") ~ 1)
ci <- data.frame(survfit.mstate.simple$prev)
names(ci) <- c("CI.1", "CI.2", "Surv")
ci


# remove when finished with these
rm(tt)
rm(status)

# from the book:
#ci.prostate <- Cuminc(time=prostateSurvival.highrisk$survTime, 
#   status=prostateSurvival.highrisk$status)

#ci1 <- ci.prostate$CI.1  # CI.1 is for prostate cancer 
#ci2 <- ci.prostate$CI.2  # CI.2 is for other causes 
#times <- ci.prostate$time/12  # convert months to years 

# use this code instead:
sf <- survfit(Surv(survTime, status, type="mstate") ~ 1, data=prostateSurvival.highrisk)
tt <- sf$time
CIs <- sf$prev
SSurv <- 1 - apply(CIs, 1, sum)
ci1 <- CIs[,1]
ci2 <- CIs[,2]
times <- tt/12

Rci2 <- 1 - ci2

# figure 9.4   (simple)       # works for R version 3.2, perhaps not for version 3.3
windows()
plot(Rci2 ~ times, type="s", ylim=c(0,1), lwd=2, col="green",
    xlab="Time in years", ylab="Survival probability",
    xlim=c(0,10))
lines(ci1 ~ times, type="s", lwd=2, col="blue") 
lines(surv.other.km ~ time.km, type="s", col="lightgreen", lwd=1)
lines(cumDist.prost.km ~ time.km, type="s", col="lightblue", lwd=1) 

# Figure 9.4 (more complete)    # works for R version 3.2, perhaps not for version 3.3
windows(height=5, width=7)

layout(matrix(c(2,1,3), ncol=3),
  heights=c(5,5,5), widths=c(0.5, 5, 0.5))
layout.show(3)
par(mar=c(5,2,2,2), cex.axis=1.65, cex.lab=1.65)

#par(mar=c(5, 4, 4, 4) + 0.1)
plot(cumDist.prost.km ~ time.km, type="s", ylim=c(0,1), lwd=1, xlab="Years from prostate cancer diagnosis",
   axes=F, col="lightblue", xlim=c(0,10))
lines(surv.other.km ~ time.km, type="s", col="lightgreen", lwd=1)
axis(1)
axis(2, at=seq(0,1,0.2), las=1)
#axis(3)
axis(4, at=seq(0,1,0.2), labels=seq(1,0,-0.2), las=1)
box()
text(20/12, 0.27, label="Death from\nprostate cancer", cex=1.65)
text(20/12, 0.73, label="Death from\nother causes", cex=1.65)

lines(Rci2 ~ times, type="s", ylim=c(0,1), lwd=2, col="green")
lines(ci1 ~ times, type="s", lwd=2, col="blue")

x <- c(0,5)
y <- c(0,0.25)
par(mar=c(0,0,0,0))
plot(x ~ y,type="n",axes=F)
text(x=0.1,y=2.5,"Probability of death from prostate cancer",srt=90,cex=1.65)
plot(x ~ y,type="n",axes=F)
text(x=0.15,y=2.5,"Probability of death from other causes",srt=-90,cex=1.65)

# Figure 9.5        # works for R version 3.2, perhaps not for version 3.3
windows(height=5, width=7)

par(mar=c(5, 4, 2, 2) + 0.1, cex.lab=1.1, cex.axis=1.1)  # default

plot(ci1 ~ times, type="s", ylim=c(0,1), lwd=2, xlab="Years from prostate cancer diagnosis",
   ylab="Probability patient has died", axes=F, col="blue", xlim=c(0,10))
ci1.ci2 <- ci1 + ci2
lines(ci1.ci2 ~ times, type="s", lwd=2, col="green")

axis(1)
axis(2, at=seq(0,1,0.2), labels=seq(0,1,0.2), las=1)
#axis(3)
#axis(4)
box()
text(8, 0.1, "Death from\nprostate cancer", cex=1.1)
text(8, 0.52, "Death from\nother causes", cex=1.1)

# # # # # # # # # # #
# Section 9.2.4 regression methods for cause-specific hazards
# # # # # # # # # # #

prostateSurvival.T2 <- prostateSurvival[prostateSurvival$stage=="T2",]
attach(prostateSurvival.T2)
result.prostate <- coxph(Surv(survTime, status.prost) ~ grade + ageGroup)
summary(result.prostate) 

result.other <- coxph(Surv(survTime, status.other) ~ grade + ageGroup)
summary(result.other) 

cov.matrix <- model.matrix(~ grade + ageGroup)
head(cov.matrix)

install.packages("cmprsk") # must be done once
library(cmprsk)
result.prostate.crr <- crr(survTime, status, cov1=cov.matrix[,-1], failcode=1) 
summary(result.prostate.crr)

result.other.crr <- crr(survTime, status, cov1=cov.matrix[,-1], failcode=2) 
summary(result.other.crr)

# # # # # # # # # # # # # 
# Section 9.2.5 covariates on different causes of death
# # # # # # # # # # # # #

tmat <- trans.comprisk(2, names = c("event-free", "prostate", "other"))
tmat

prostate.long <- msprep(time = cbind(NA, survTime, survTime), 
      status = cbind(NA, status.prost, status.other),  
      keep = data.frame(grade, ageGroup), trans = tmat) 
head(prostate.long) 
events(prostate.long) 

summary(coxph(Surv(time, status) ~ grade + ageGroup, data=prostate.long, subset={trans==1}))
summary(coxph(Surv(time, status) ~ grade + ageGroup, data=prostate.long, subset={trans==2}))

summary(coxph(Surv(time, status) ~ grade + ageGroup + strata(trans), data=prostate.long))
summary(coxph(Surv(time, status) ~ grade*factor(trans) + ageGroup + strata(trans), data=prostate.long))  
summary(coxph(Surv(time, status) ~ (grade + ageGroup)*trans + ageGroup + strata(trans), data=prostate.long)) 

