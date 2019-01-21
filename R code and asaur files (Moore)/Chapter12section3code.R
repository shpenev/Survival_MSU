# # # # # # # # # # # # # # # # # # #
# Applied Survival Analysis Using R
# Dirk F. Moore
# Springer (2016)
# # # # # # # # # # # # # # # # # # #

# Chap. 12 Sect. 12.3  lasso for predictive markers

# install.packages("asaur")   # must be done once
install.packages("penalized") # must be done once

library(asaur)
hepatoCellularNoMissing <- 
    hepatoCellular[complete.cases(hepatoCellular),] 
hepatoCellularNoMissing[c(1,5,12),c(16,17, 23:27)]

attach(hepatoCellularNoMissing)
library(penalized)
hepato.pen <- penalized(Surv(OS, Death),  
    penalized=hepatoCellularNoMissing[,23:48],  
    standardize=T, lambda1=10) 
round(coef(hepato.pen, standardize=T), 3)

set.seed(34)
hepato.prof <- profL1(Surv(OS, Death),
   penalized=hepatoCellularNoMissing[,23:48],
   standardize=T, fold=10, minlambda1=2, maxlambda1=12)

# Fig. 12.7 partial likelihood
plot(hepato.prof$cvl ~ hepato.prof$lambda, type="l", log="x",
     xlab="lambda", ylab="Cross-validated log partial likelihood")

set.seed(34) 
hepato.opt <- optL1(Surv(OS, Death),    
 	penalized=hepatoCellularNoMissing[,23:48],  standardize=T, fold=10) 
hepato.opt$lambda 

# Fig. 12.7  add vertical line at maximum
abline(v=hepato.opt$lambda, col="gray")

hepato.pen <- penalized(Surv(OS, Death), 
  penalized=hepatoCellularNoMissing[,23:48],  standardize=T, 
  steps=20, lambda1=8)

# Fig. 12.8
plotpath(hepato.pen, labelsize=0.9, standardize=T, log="x", lwd=2) 
abline(v=hepato.opt$lambda, col="gray", lwd=2)

hepato.pen <- penalized(Surv(OS, Death),
    penalized=hepatoCellularNoMissing[,23:48],  standardize=T,
    lambda1=hepato.opt$lambda) 
round(coef(hepato.pen, standardize=T), 3)

# Figure 12.9  without labels
plot(predict(hepato.pen, 
   penalized=hepatoCellularNoMissing[c(1, 5, 12),23:48]))

hepato.predict.1 <- predict(hepato.pen,
   penalized=hepatoCellularNoMissing[1,23:48]) 
hepato.predict.5 <- predict(hepato.pen,
   penalized=hepatoCellularNoMissing[5,23:48]) 
hepato.predict.12 <- predict(hepato.pen, 
   penalized=hepatoCellularNoMissing[12,23:48])

slotNames(hepato.predict.1)
# Figure 12.9

par(mar=c(5, 5, 4, 2) + 0.1)
plot(stepfun(hepato.predict.1@time[-1], hepato.predict.1@curves), do.points=F,
  ylim=c(0,1), xlab="Time in months", ylab="Predicted survival probability",
  cex.axis=1.4, cex.lab=1.4, lwd=2, main="", xlim=c(0,50))
plot(stepfun(hepato.predict.5@time[-1], hepato.predict.5@curves), 
  do.points=F, add=T, col="blue", lwd=2)
plot(stepfun(hepato.predict.12@time[-1], hepato.predict.12@curves), do.points=F, add=T, col="red", lwd=2)

legend("bottomleft", legend=c("Patient 1", "Patient 5", "Patient 12"), 
  col=c("black", "blue", "red"), lwd=2, cex=1.4)

 