# # # # # # # # # # # # # # # # # # #
# Applied Survival Analysis Using R
# Dirk F. Moore
# Springer (2016)
# # # # # # # # # # # # # # # # # # #

# Write out the data sets in the "asaur" package as comma-separated text files
#  for input to other statistical packages

# Source:  "Applied Survival Analysis Using R" by Dirk F. Moore  Springer, 2016
install.packages("asaur")   # this must be done once
setwd("C:\\asaur") # point to a pre-defined directory where the data will be written; 
                   # note use of double "\\" to define a Windows directory
                   # Use "/" for Macintosh or unix operating system
library(asaur)

write.csv(ashkenazi, file="ashkenazi.csv", row.names=F)

write.csv(ChanningHouse, file="ChanningHouse.csv", row.names=F)

write.csv(gastricXelox, file="gastricXelox.csv", row.names=F)

write.csv(hepatoCellular, file="hepatoCellular.csv", row.names=F, na=".")

write.csv(pancreatic, file="pancreatic.csv", row.names=F)

write.csv(pancreatic2, file="pancreatic2.csv", row.names=F)

write.csv(pharmacoSmoking, file="pharmacoSmoking.csv", row.names=F)

write.csv(prostateSurvival, file="prostateSurvival.csv", row.names=F)
