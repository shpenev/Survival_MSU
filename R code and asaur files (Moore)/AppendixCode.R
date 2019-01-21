# # # # # # # # # # # # # # # # # # #
# Applied Survival Analysis Using R
# Dirk F. Moore
# Springer (2016)
# Appendix
# # # # # # # # # # # # # # # # # # #

# # # # # # # # # # # # #
# A.1.1 First R session
# # # # # # # # # # # # #

x <- 3
y <- 2
y^x
y**x

x.vec <- c(1, 3.5, 7)
y.vec <- c(2, 7, 8.6)
x.vec
y.vec
x.vec + y.vec
z.vec <- x.vec + y.vec
z.vec

z.vec <- c(x.vec, y.vec)
z.vec
z.vec[1:4]

w.vec <- c("a", "A", "aBc")
w.vec

library(survival)
tt <- c(2, 5, 6, 7, 8)
status <- c(1, 1, 0, 1, 1)
Surv(tt, status)

tt <- c(2, 5, 6,
   7, 8)
tt

Surv(tt, status) # Create a survival data structure

# # # # # # # # # # # #
# quit R session if you like
q()

# # # # # # # # # # # # #
# A.1.2 Scatterplots and linear regression
# # # # # # # # # # # # #

x.vec <- 1:10 
x.vec 

y.vec <- 3 + 2*x.vec + rnorm(10, mean=0, sd=2) 
y.vec   
plot(y.vec ~ x.vec)

plot(y.vec ~ x.vec, xlim=c(0, 10), ylim=c(0, 30),
    xlab="x", ylab="y") 
title("A simple plot")

result.lm <- lm(y.vec ~ x.vec)
result.lm
abline(result.lm)

# # # # # # # # # # # # #
# A.1.3 Non-linear relationships
# # # # # # # # # # # # #

ff <- function(x) {
  result <- 2*x^3 - 9*x^2 + 5*x + 6
  result   
  }
ff(x=c(0, 1, 2))

x.vec <- (-99:400)/100    # create 500 points between -1 and 4
y.vec <- ff(x.vec) + 10*rnorm(500)   # fixed and random effects

# Fig. A.1
plot(y.vec ~ x.vec, col="gray")
curve(ff, from=-1, to=4, col="red", lwd=2, add=T) 

x2.vec <- x.vec^2
x3.vec <- x.vec^3
result.lm <- lm(y.vec ~ x.vec + x2.vec + x3.vec)
summary(result.lm)

result.smooth <- loess(y.vec ~ x.vec) 
smooth.estimates <- predict(result.smooth) 
lines(smooth.estimates ~ x.vec, col="blue", lwd=2)

smoothSEcurve <- function(yy, xx) {
  # use after a call to "plot"   
  # fit a lowess curve and 95% confidence interval curve
  
  # make list of x values
  xx.list <- min(xx) + ((0:100)/100)*(max(xx) - min(xx)) 

  # Then fit loess function through the points (xx, yy)
  #   at the listed values
  yy.xx <- predict(loess(yy ~ xx), se=T, 
     newdata=data.frame(xx=xx.list))

  lines(yy.xx$fit ~ xx.list, lwd=2)
  lines(yy.xx$fit - 
       qt(0.975, yy.xx$df)*yy.xx$se.fit ~ xx.list, lty=2)
  lines(yy.xx$fit + 
       qt(0.975, yy.xx$df)*yy.xx$se.fit ~ xx.list, lty=2)
  }
smoothSEcurve(y.vec, x.vec)

# # # # # # # # # # # # #
# A.1.4 Data frames and search path
# # # # # # # # # # # # #

lung[1:6,1:7]   
time.A <- lung[,2]
time.B <- lung$time 
time.A[1:5] 
time.B[1:5] 

attach(lung)
time[1:5] 

time <- c(1,2,3,4) 
time
rm(time)
time[1:5] 

# # # # # # # # # # # # #
# A.1.5 Defining variables within a data frame
# # # # # # # # # # # # #

lung <- within(lung, {
     delta <- status - 1 })
lung[1:6, c(1:7, 11)]

# # # # # # # # # # # # # # # #
# Sec. A.1.6 Importing and exporting data frames
# # # # # # # # # # # # # # # #

# the following only works of the directory "c:\survival" already exists
setwd("c:\\survival")  # you must create a directory "c:\survival" if it doesn't already exist
write.csv(lung, file="lung.csv", na=".", row.names=F)

setwd("c:\\survival")
lung2 <- read.csv("lung.csv", na.strings=".", header=T)
head(lung2)

# # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # #

# # # # # # # # # # # # # # # #
# Sec. A.2.2 Working with dates in R   Using the "as.date" function
# # # # # # # # # # # # # # # #

install.packages("date")
library("date")
date.1 <- as.date("8/31/1956") 
date.2 <- as.date("7/5/1957") 
date.1 
date.2

date.2 - date.1 
as.numeric(date.1) 
as.numeric(as.date("1/1/1960")) 
as.date("2/29/2000") 
as.date("2/29/1900") 
as.date("January 30 2005") 
as.date("30/1/2005", order="dmy") 

entry.dates <- c("9/20/2010", "9/30/2010", "11/2/2010", "1/5/2011") 
death.dates <- c("5/4/2013", NA, "6/9/2013", "4/5/2012") 
lastSeen.dates <- c("5/4/2013", "8/21/2013", "6/9/2013", "4/5/2012") 
  
entry <- as.date(entry.dates) 
death <- as.date(death.dates) 
lastSeen <- as.date(lastSeen.dates) 

censor <- as.numeric(!is.na(death))
censor 

survTime.temp <- death - entry
survTime.temp 

survTime <- survTime.temp
survTime[censor == 0] <- lastSeen[censor == 0] - entry[censor == 0]
survTime 
censor 

library(survival)
Surv(survTime, censor) 

# # # # # # # # # # # # # # # #
# Sec. A.3  forest plots
# # # # # # # # # # # # # # # #

library(survival)
head(veteran)
attach(veteran)  # this is necessary to make the variables in "veteran" visible

trt.f <- factor(trt, labels=c("standard", "test"))
result <- coxph(Surv(time, status) ~ trt.f + celltype, 
    data=veteran)
result 

coef.est <- c(NA, NA, 0, 0.198, NA, NA, NA, 0, 1.096, 1.169, 0.297) 
se.est <- c(NA, NA, 0, 0.197, NA, NA, NA, 0, 0.272, 0.295, 0.286) 
lower <- coef.est - 1.96*se.est 
upper <- coef.est + 1.96*se.est
label.factors <- matrix(c("Treatment Group", "", "  standard",
   "  test", "", "Cell Type", "", "  sqamous", "  smallcell", 
   "  adeno", "  large"), ncol=1)

# Fig. A.2

install.packages("forestplot")    # this must be done once
library(forestplot) 
windows()
forestplot(label.factors, coef.est,  lower=lower, upper=upper,
      boxsize=0.4, xticks=c(-0.5,0,0.5, 1, 1.5, 2),
      txt_gp=fpTxtGp(label=gpar(cex=1.5)))

# # # # # # # # # # # # # # # #
# Sec. A.4  extracting the log partial likelihood and estimates from a coxph object
# # # # # # # # # # # # # # # #

rm(status)   # delete this old variable

library(survival)
attach(veteran) 
testInd <- trt - 1  # now 0 refers to standard, and 1 to test
result <- coxph(Surv(time, status) ~ testInd)
result 

result.cox.0 <- coxph(Surv(time, status) ~ testInd, 
     init=0, control=list(iter.max=0))
loglik.0 <- result.cox.0$loglik[2]
loglik.0

coef.mple <- as.numeric(result$coef)
coef.mple 

result.cox.max <- coxph(Surv(time, status) ~ testInd, 
     init=coef.mple, control=list(iter.max=0))
loglik.max <- result.cox.max$loglik[2]
2*(loglik.max - loglik.0) 

pchisq(0.009643379, 1, lower.tail=F) 

 


