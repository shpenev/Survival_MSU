# # # # # # # # # # # # # # # # # # #
# Applied Survival Analysis Using R
# Dirk F. Moore
# Springer (2016)
# # # # # # # # # # # # # # # # # # #

# Chapter 4

# install.packages("asaur")   # must do this once
library(survival)

# # # # # # # # # #
# Section 4.1, simple example
# # # # # # # # # #

tt <- c(6, 7, 10, 15, 19, 25)
delta <- c(1, 0, 1, 1, 0, 1)
trt <- c(0, 0, 1, 0, 1, 1)

survdiff(Surv(tt, delta) ~ trt)

# # # # # # # # # #
# Section 4.1, pancreatic data
# # # # # # # # # #

install.packages("date")  # must be done once
library(date)     # load the "date" library
library(survival) # load the "survival" library
library(asaur)

head(pancreatic)

attach(pancreatic)     # make the variable names accessible

# convert the text dates into R dates
Progression.d <- as.date(as.character(progression))
OnStudy.d <- as.date(as.character(onstudy))
Death.d <- as.date(as.character(death))

# compute progression-free survival

progressionOnly <- Progression.d - OnStudy.d
overallSurvival <- Death.d - OnStudy.d
pfs <- pmin(progressionOnly, overallSurvival)
pfs[is.na(pfs)] <- overallSurvival[is.na(pfs)]

# convert pfs to months

os <- overallSurvival
status <- rep(1, length(pfs))

# Note:  "pancreatic2" is already part of the "asaur" package; this is how one would create it:
#pancreatic2 <- data.frame(pfs, os, status, stage)    # create "pancreatic2" data frame

# library(pancreatic2) # an alternative to the above is to just use pancreatic2
pfs.month <- pfs/30.5

plot(survfit(Surv(pfs.month) ~ stage), xlab="Time in months",
   ylab="Survival probability", col=c("blue", "red"), lwd=2)
legend("topright", legend=c("Locally advanced", "Metastatic"),
  col=c("blue","red"), lwd=2)
  

result.survdiff.0 <- survdiff(Surv(pfs) ~ stage, rho=0)
result.survdiff.0    # not statistically significant

# Carry out the Prentice modification of the Gehan test (rho=1)
result.survdiff.1 <- survdiff(Surv(pfs) ~ stage, rho=1)
result.survdiff.1  #statistically significant

# # # # # # # # # # # #
# pharmacoSmoking example, Section 4.2
# # # # # # # # # # # #

library(asaur)

attach(pharmacoSmoking)
library(survival)

survdiff(Surv(ttr, relapse) ~ grp)

table(ageGroup2)
survdiff(Surv(ttr, relapse) ~ grp + strata(ageGroup2))

# # # # # # # # # # # #
# genotype confounding simulation, Section 4.2, Example 4.3
# # # # # # # # # # # #

lambda.mutant.0 <- 0.03 
lambda.mutant.1 <- 0.03*0.55 
lambda.wt.0 <- 0.03*0.2 
lambda.wt.1 <- 0.03*0.2*0.55

set.seed(4321)

tt.control.mutant <- rexp(25, rate=lambda.mutant.0) 
tt.treat.mutant <- rexp(125, rate=lambda.mutant.1) 
tt.control.wt <- rexp(125, rate=lambda.wt.0) 
tt.treat.wt <- rexp(25, rate=lambda.wt.1)
ttAll <- c(tt.control.mutant, tt.treat.mutant, tt.control.wt, 
    tt.treat.wt)

status <- rep(1, length(ttAll)) 

genotype <- c(rep("mutant", 150), rep("wt", 150)) 
trt <- c(rep(0, 25), rep(1, 125), rep(0, 125), rep(1, 25))

survdiff(Surv(ttAll, status) ~ trt) 
survdiff(Surv(ttAll, status) ~ trt + strata(genotype))
geneConfounder <- data.frame(ttAll, status, trt, genotype)

# # # # # # # # # #
# Section 4.2, Figure 4.2.1
# # # # # # # # # #

windows()

plot(survfit(Surv(ttAll, status) ~ trt, data=geneConfounder),
  col=c("red", "blue"), xlab="Time in days", xlim=c(0,800),
  ylab="Survival probability", lwd=2, cex.axis=1.2, cex.lab=1.2)
legend("topright", 
    legend=c("control", "treated"),
    lty=c(1, 1), col=c("red", "blue"), lwd=2, cex=1.2)
    
windows()    
plot(survfit(Surv(ttAll, status) ~ trt, data=geneConfounder,
  subset={genotype == "mutant"}), xlim=c(0,800),
  col=c("red", "blue"), xlab="Time in days",
  ylab="Survival probability", lwd=2, cex.axis=1.2, cex.lab=1.2)

lines(survfit(Surv(ttAll, status) ~ trt, data=geneConfounder, 
   subset={genotype == "wt"}), 
   lty=2,  lwd=2, col=c("red", "blue"))

legend("topright", 
    legend=c("control mutant", "treated mutant", "control wildtype", "treated wildtype"),
    lty=c(1, 1, 2, 2), col=c("red", "blue", "red", "blue"), lwd=2, cex=1.2)



