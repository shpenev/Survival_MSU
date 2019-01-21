# # # # # # # # # # # # # # # # # # #
# Applied Survival Analysis Using R
# Dirk F. Moore
# Springer (2016)
# # # # # # # # # # # # # # # # # # #

# Chapter 6

# install.packages("asaur")    # must be done once

# # # # # # # # # # # #
# repeat generation of genotype simulation
# genotype confounding simulation, Section 4.2, Example 4.3
# # # # # # # # # # # #

library(survival)
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

geneConfounder <- data.frame(ttAll, status, trt, genotype)

# # # # # # # # # # # #
# Now use this simulated data to fit model in Section 6.1
# # # # # # # # # # # #

coxph(Surv(ttAll, status) ~ trt)

coxph(Surv(ttAll, status) ~ trt + strata(genotype))

coxph(Surv(ttAll, status) ~ trt + genotype)

# # # # # # # # # # # #
# Section 6.2   Categorical and continuous covariates
# # # # # # # # # # # #

race <- factor(c("black", "black", "white", "white", "other", 
           "other"))
age <- c(48, 52, 87, 82, 67, 53)
model.matrix(~ race + age)[,-1]

race <- relevel(race, ref="white")
model.matrix(~ race + age)[,-1]
model.matrix(~ race + age + race:age)[,-1]   

age <- runif(n=60, min=40, max=80) 
race <- factor(c(rep("white", 20), rep("black", 20), 
          rep("other", 20))) 
race <- relevel(race, ref="white") 
log.rate.vec <- -4.5 + c(rep(0,20), rep(1,20), rep(2,20)) + age*0.05 
tt <- rexp(n=60, rate=exp(log.rate.vec)) 
status <- rep(1, 60)
library(survival)
result.cox <- coxph(Surv(tt, status) ~ race + age) 
summary(result.cox)

# remove these created variables when finished
rm(race)
rm(age)
rm(tt)
rm(status)
# # # # # # # # # # # #
# Section 6.3    hypothesis testing for nested models
# # # # # # # # # # # #

library(asaur)
library(survival)
attach(pharmacoSmoking)
levels(ageGroup4)
levels(employment) 

modelA.coxph <- coxph(Surv(ttr, relapse) ~ ageGroup4)
modelA.coxph

modelB.coxph <- coxph(Surv(ttr, relapse) ~ employment)
modelB.coxph

modelC.coxph <- coxph(Surv(ttr, relapse) ~ ageGroup4 + employment)
modelC.coxph 

logLik(modelA.coxph)
logLik(modelB.coxph)
logLik(modelC.coxph)

pchisq(4.567, df=2, lower.tail=F) 
pchisq(14.727, df=3, lower.tail=F) 

model.null.coxph <- coxph(Surv(ttr, relapse) ~ 1)
logLik(model.null.coxph) 
pchisq(12.2206, df=3, lower.tail=F)

anova(modelA.coxph, modelC.coxph) 

# forest plots, Fig 6.1

coef.est <- c(NA, 0, 0.656, NA, NA, 0, 0.623, 0.521,
           NA, NA, 0, -0.112, -1.023, -0.707)
se.est <- c(NA, 0, 0.220, NA, NA, 0, 0.276, 0.332,
           NA, NA, 0, 0.332, 0.360, 0.502)
lower <- coef.est - 1.96*se.est
upper <- coef.est + 1.96*se.est


label.factors <- c("Treatment Group", "   triple therapy", "   patch", "",
   "Employment", "   full time", "   other", "   part time",
   "", "Age group", "   21-34", "   35-49", "   50-64", "   65+")
   
install.packages("forestplot")  # must do this once
library(forestplot)

forestplot(label.factors,
  coef.est,
  lower,
  upper,
  zero = 0,
  cex = 1.0,
  lineheight = "auto",
  xlab = "Log hazard ratio",
  boxsize=0.4,
  xticks=c(-1.5, -1.0, -0.5,0,0.5, 1, 1.5),
  txt_gp=fpTxtGp(label=gpar(cex=1.3)),
  new_page=T)

# # # # # # # # # # # #
# Section 6.4   Akaike Information Criterion for non-nested models
# # # # # # # # # # # #

AIC(modelA.coxph) 
AIC(modelB.coxph) 
AIC(modelC.coxph) 

# be sure you have deleted "race" and "age" before running the following!
modelAll.coxph <- coxph(Surv(ttr, relapse) ~ grp + gender + race +
            employment + yearsSmoking + levelSmoking + ageGroup4 +
            priorAttempts + longestNoSmoke)

result.step <- step(modelAll.coxph, scope=list(upper=~ grp + 
            gender + race + employment + yearsSmoking + 
            levelSmoking + ageGroup4 + priorAttempts + 
            longestNoSmoke, lower=~grp) )
result.step

# # # # # # # # # # # #
# Section 6.5  
# # # # # # # # # # # #
            
modelS4.coxph <- coxph(Surv(ttr, relapse) ~ grp + employment + pspline(age, df=4) ) 
modelS4.coxph
 
termplot(modelS4.coxph, se=T, terms=3, ylabs="Log hazard") 
 