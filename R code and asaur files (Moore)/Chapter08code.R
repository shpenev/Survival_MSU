# # # # # # # # # # # # # # # # # # #
# Applied Survival Analysis Using R
# Dirk F. Moore
# Springer (2016)
# # # # # # # # # # # # # # # # # # #

# Chapter 8

# install.packages("asaur")    # this must be done once

# Section 8.1 Stanford heart transplant data
library(survival)
coxph(Surv(futime, fustat) ~ transplant + age + surgery, data=jasa) # naive analysis

# landmark method, with landmark set to 30 days
ind30 <- jasa$futime >= 30
transplant30 <- {{jasa$transplant == 1} & {jasa$wait.time < 30}}
coxph(Surv(futime, fustat) ~ transplant30 + age + surgery, data=jasa, subset=ind30)

# # # # # # # # # # #
# subset example
# # # # # # # # # # #

id <- 1:nrow(jasa)
jasaT <- data.frame(id, jasa)
id.simple <- c(2, 5, 10, 12, 28, 95)
# jasaT[jasaT$futime < 100,c(6,8,9,10)]
# heart.simple <- jasaT[c(2, 5, 10, 12, 74, 95),c(1, 10, 9, 6, 11)]
heart.simple <- jasaT[id.simple,c(1, 10, 9, 6, 11)]

coxph(Surv(futime, fustat) ~ transplant, data=heart.simple) # naive analysis

# put data into counting process format
sdata <- tmerge(heart.simple, heart.simple, id=id,
                death=event(futime, fustat),
                transpl=tdc(wait.time))
heart.simple.counting <- sdata[,-(2:5)]     
heart.simple.counting

coxph.heart.simple.counting <- coxph(Surv(tstart, tstop, death) ~ transpl,
   data=heart.simple.counting) 
summary(coxph.heart.simple.counting)   # time-dependent analysis 

# # # # # # # # # # #
# counting process format for entire Stanford heart transplant data
# See Therneau and Crowson (2015) Using time dependent covariates and 
#   time-dependent coefficients in the Cox model. 
#   Vignette for R survival package
#   http://cran.r-project.org/web/packages/survival
# # # # # # # # # # #

tdata <- jasa[, -c(1:4, 11:14)]  #leave off the dates and transplant-specific covariates, temporary data set
tdata$futime <- pmax(.5, tdata$futime)  # the death on day 0
indx <- {{tdata$wait.time == tdata$futime} & !is.na(tdata$wait.time)}
#indx <- with(tdata, which(wait.time == futime))
tdata$wait.time[indx] <- tdata$wait.time[indx] - .5  #the tied transplant
id <- 1:nrow(tdata)
tdata$id <- id
sdata <- tmerge(tdata, tdata, id=id, 
                death = event(futime, fustat), 
                trans = tdc(wait.time))
jasa.counting <- sdata[,c(7:11, 2:3)]
head(jasa.counting)
summary(coxph(Surv(tstart, tstop, death) ~ trans + surgery + age, data=jasa.counting))

# # # # # # # # # # # #
# Section 8.2.1 predictable time-dependent variables
# # # # # # # # # # # #

library(asaur)

attach(pancreatic2)

stage.n <- rep(0, nrow(pancreatic2))
stage.n[pancreatic2$stage == "M"] <- 1

result.panc <- coxph(Surv(pfs) ~ stage.n)   # this is the log-rank test
result.panc

result.panc2.tt <- coxph(Surv(pfs) ~ stage.n + tt(stage.n),
    tt=function(x,t, ...) x*log(t))
result.panc2.tt

result.sch.resid <- cox.zph(result.panc, transform=function(pfs) log(pfs))
plot(result.sch.resid, cex.axis=1.5, cex.lab=1.5)
abline(coef(result.panc2.tt), col="red")

# # # # # # # # # # # #
# Section 8.2.2 
# # # # # # # # # # # #

coxph(Surv(time, status==2) ~ age, data=lung)

# the following enters age as time-dependent, but gives the same result:
coxph(Surv(time, status==2) ~ tt(age), data=lung,
  tt=function(x, t, ...) {
    age <- x + t/365.25
    age})
    
    

