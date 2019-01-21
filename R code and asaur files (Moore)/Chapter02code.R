# # # # # # # # # # # # # # # # # # #
# Applied Survival Analysis Using R
# Dirk F. Moore
# Springer (2016)
# # # # # # # # # # # # # # # # # # #

# code for Chapter 2

library(survival)
tm <- c(0, # birth
   1/365, # first day of life
   7/365, # seventh day of life
   28/365, # fourth week of life
   1:110) # subsequent years
hazMale <- as.numeric(survexp.us[,"male","2004"]) # 2004 males
hazFemale <- as.numeric(survexp.us[,"female","2004"]) # 2004 females
tm.diff <- diff(tm)
survMale <- exp(-cumsum(hazMale*tm.diff)*365.24)
survFemale <- exp(-cumsum(hazFemale*tm.diff)*365.24)

# Figure 2.1.2 log hazard and survival for US males and females in 2004
par(mfrow=c(2,1),    # two rows and one column of plots
    mar=c(4.2,5,2,2))  # set margins for the lower, left, top, and righ of each plot

logHazMale <- log(hazMale)
logHazFemale <- log(hazFemale)    

plot(logHazMale ~ tm[-1], type="l",
     xlab="Age in years",           # x axis label
     ylab="Hazard",col="blue",      # y azis label
     lwd=2,                         # double line width
     las=1,                         # make y axis labels perpendicular to axis
     axes=F, cex.lab=1.3, cex.axis=1.3)     # make blue line solid
lines(logHazFemale ~ tm[-1],type="l", 
      col="red",lwd=2, lty=2)   # add a red dashed line to the plot

yyLabs <- c(1e-07, 1e-06, 1e-05, 1e-04, 1e-03, 1e-02)
yyLabsLog <- log(yyLabs)
axis(2, at=yyLabsLog, labels=c(expression(10^-7), expression(10^-6), 
  expression(10^-5), expression(10^-4), expression(10^-3), expression(10^-2)), las=1)  
axis(1, cex.axis=1.3)   
legend("bottomright", legend=c("males","females"),
       lty=c(1,2), col=c("blue","red"), lwd=2, cex=1.3)
title("Hazards for US males and females in 2004")

tm.diff <- diff(tm)         # same length as "tm"
survMale <- exp(-cumsum(hazMale*tm.diff)*365.24)         # survival probs for males
survFemale <- exp(-cumsum(hazFemale*tm.diff)*365.24)     # survival probs for females
#windows(width=7,height=5)
plot(survMale ~ tm[-1],type="l",          # lower case "L" indicates line plot
     xlab="Age in years",             # x axis label
     ylab="Survival probability",     # y azis label
     col="blue",                      # line color
     lwd=2,                           # double line width
     las=1,                           # make y axis labels perpendicular to axis
     ylim=c(0,1), cex.lab=1.3, cex.axis=1.3)       # y axis limit ranges from 0 to 1

lines(survFemale ~ tm[-1], col="red", lwd=2, lty=2)    # add a red dashed line to the plot
legend("bottomleft", legend=c("males","females"),
       lty=c(1,2), col=c("blue","red"), lwd=2, cex=1.3)
title("Survival of US males and females in 2004")


# mean age of death for men and women
sum(survMale*tm.diff) # mean age of male death in 2004
sum(survFemale*tm.diff) # mean age of female death in 2004

# # # # #
# Figure 2.4.2 Weibull hazard functions
# # # # #

windows()
weibHaz <- {function(x, shape, scale) dweibull(x, shape=shape,
     scale=scale)/pweibull(x, shape=shape, scale=scale, lower.tail=F)}
curve(weibHaz(x, shape=1.5, scale=1/0.03), from=0, to=80, 
     ylab='Hazard', xlab='Time', col="blue")
curve(weibHaz(x, shape=1.0, scale=1/0.03), from=0, to=80, 
     ylab='Hazard', xlab='Time', add=T, col="black")
curve(weibHaz(x, shape=0.75, scale=1/0.03), from=0, to=80, 
     ylab='Hazard', xlab='Time', add=T, col="red")
text(45, 0.065, expression(alpha == 1.50), col="red", cex=1.3)
text(58, 0.065, expression(lambda==0.03), col="red", cex=1.3)
text(45, 0.015, expression(alpha == 0.75), col="blue", cex=1.3)
text(58, 0.015, expression(lambda==0.03), col="blue", cex=1.3)
text(45, 0.034, expression(alpha == 1.00), col="black", cex=1.3)
text(58, 0.034, expression(lambda==0.03), col="black", cex=1.3)

# # # # #
# Figure 2.4.3 Gamma hazard functions
# # # # #

windows()
gammaHaz <- {function(x, shape, scale) dgamma(x, shape=shape,
     scale=scale)/pgamma(x, shape=shape, scale=scale, lower.tail=F)}
curve(gammaHaz(x, shape=1.5, scale=1/0.03), from=0, to=80, 
     ylab='Hazard', xlab='Time', col="red", ylim=c(0,0.08) )
curve(gammaHaz(x, shape=0.75, scale=1/0.03), from=0, to=80, 
     ylab='Hazard', xlab='Time', col="blue", add=T )
curve(gammaHaz(x, shape=1.0, scale=1/0.03), from=0, to=80, 
     ylab='Hazard', xlab='Time', col="black", add=T )
text(12, 0.01, expression(beta == 1.50), col="red", cex=1.3)
text(25, 0.01, expression(lambda==0.03), col="red", cex=1.3)
text(12, 0.05, expression(beta == 0.75), col="blue", cex=1.3)
text(25, 0.05, expression(lambda==0.03), col="blue", cex=1.3)
text(12, 0.027, expression(beta == 1.00), col="black", cex=1.3)
text(25, 0.027, expression(lambda==0.03), col="black", cex=1.3)
