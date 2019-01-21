# # # # # # # # # # # # # # # # # # #
# Applied Survival Analysis Using R
# Dirk F. Moore
# Springer (2016)
# # # # # # # # # # # # # # # # # # #

# Chapter 3 code

# Kaplan-Meier estimate of the survival function using the data
#   in Example 1.1

# install.packages("asaur")  # must do this once
library(survival)
tt <- c(7,6,6,2,4)
cens <- c(0,1,0,1,1)
Surv(tt, cens)
result.km <- survfit(Surv(tt, cens) ~ 1, conf.type="log-log")
summary(result.km)

# Figure 3.1.1 Kaplan-Meier plot
plot(result.km, conf.int=F, ylab="Survival probability",
      xlab="Time in years")

# Figure 3.1.1 nicer plot
plot(result.km, conf.int=F, ylab="Survival probability",
      xlab="Time in years", cex.lab=1.5, cex.axis=1.5, lwd=2)

# Figure 3.1.2 Kaplan-Meier survival curve with 95% confidence intervals
plot(result.km, conf.int=T, ylab="Survival probability",
      xlab="Time in years")

# # # # #      
# Figure 3.1.3   Kaplan-Meier progression-free survival plot for gastricXelox data
# # # # #
library(asaur)
library(survival)

timeMonths <- gastricXelox$timeWeeks*7/30.25
result.km <- survfit(Surv(timeMonths, delta) ~ 1, conf.type="log-log", data=gastricXelox)
plot(result.km, mark="|", ylab="Survival probability", xlab="Time in months",
     cex.axis=1.5, cex.lab=1.5, lwd=1.5)

# median survival and 95% confidence interval is printed as follows:
result.km

# median follow-up time
delta.followup <- 1 - gastricXelox$delta
survfit(Surv(timeMonths, delta.followup) ~ 1)

# # # # # #
# Figure 3.4.1 smooth hazard estimate
# # # # # #

install.packages("muhaz")    # this must be done once
library(muhaz)
t.vec <- c(7,6,6,5,2,4)
cens.vec <- c(0,1,0,0,1,1)

result.simple <- muhaz(t.vec, cens.vec, max.time=8,
     bw.grid=2.25, bw.method="global", b.cor="none")
plot(result.simple, xlab="Time", ylab="Hazard", lwd=2, col="blue")

# # # # #
# Figure 3.4.2 step versus smooth hazard estimate
# # # # #



result.pe5 <- pehaz(timeMonths, gastricXelox$delta, width=5, max.time=20)
plot(result.pe5, ylim=c(0,0.15), col="black")
result.pe1 <- pehaz(timeMonths, gastricXelox$delta, width=1, max.time=20)
lines(result.pe1)
result.smooth <- muhaz(timeMonths, gastricXelox$delta, bw.smooth=20,
     b.cor="left", max.time=20)
lines(result.smooth)

# # # # # #
# Figure 3.4.3 Kaplan-Meier survival plot and smoothed survival estimate for gastricXelox data
# # # # # #

haz <- result.smooth$haz.est
times <- result.smooth$est.grid
surv <- exp(-cumsum(haz[1:(length(haz)-1)]*diff(times)))
result.km <- survfit(Surv(timeMonths, gastricXelox$delta) ~ 1,
conf.type="none")
plot(result.km, conf.int=T, mark="|", xlab="Time in months",
xlim=c(0,30), ylab="Survival probability")
lines(surv ~ times[1:(length(times) - 1)])

# # # # #
# Left truncation simple example, Section 3.5
# # # # #

tt <- c(7, 6, 6, 5, 2, 4)
status <- c(0, 1, 0, 0, 1, 1)
backTime <- c(-2, -5, -3, -3, -2, -5)

tm.enter <- -backTime
tm.exit <- tt - backTime

result.left.trunc.km <- survfit(Surv(tm.enter, tm.exit, status, 
   type="counting") ~ 1, conf.type="none")
summary(result.left.trunc.km)

result.left.trunc.naa <- survfit(Surv(tm.enter, tm.exit, status, 
   type="counting") ~ 1, type="fleming-harrington", conf.type="none")
summary(result.left.trunc.naa)

# # # # # # #
# Channing House example, Section 3.5
# # # # # # #

library(asaur)

# define entry and exit times
ChanningHouse <- within(ChanningHouse, {
  entryYears <- entry/12
  exitYears <- exit/12})
head(ChanningHouse)  
ChanningMales <- ChanningHouse[ChanningHouse$sex == "Male",]
result.km <- survfit(Surv(entryYears, exitYears, cens, type="counting") ~ 1, data=ChanningMales)

# # # # # #
# Figure 3.5.3 Channing House males survival estimates
# # # # # #

plot(result.km, xlim=c(64, 101), xlab="Age", ylab="Survival probability",
  conf.int=F,
  cex.axis=1.3, cex.lab=1.3)
result.naa <- survfit(Surv(entryYears, exitYears, cens, type="counting") ~ 1, type="fleming-harrington", data=ChanningMales)
lines(result.naa, col="blue", conf.int=F)
result.km.68 <- survfit(Surv(entryYears, exitYears, cens, type="counting") ~ 1, start.time=68, data=ChanningMales)
lines(result.km.68, col="green", conf.int=F)
legend("topright", legend=c("KM", "NAA", "KM 68 and older"), lty=1, col=c("black", "blue", "green"))








