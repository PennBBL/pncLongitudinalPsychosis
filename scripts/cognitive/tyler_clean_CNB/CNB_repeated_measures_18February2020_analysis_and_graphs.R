
library(ggplot2)
library(gridExtra)
library(gtools)
library(mgcv)
library(visreg)
library(lubridate)
library(psych)

# graph main age effects

x <- read.csv("CNB_Longitudinal_Core_11February2020.csv")

x <- x[which(x$Age < 40),]

gam1 <- gamm(ACC_raw~Test+s(Age,by=Test),random=list(bblid=~1),data=x)
gam2 <- gamm(RT_raw~Test+s(Age,by=Test),random=list(bblid=~1),data=x)

gam1$gam$data <- x
gam2$gam$data <- x

pdf("CNB_repeated-measures_Age_Effects_by_test.pdf",width=12,height=6)
visreg(gam1$gam,xvar="Age",by="Test",overlay=TRUE,ylab="Accuracy",partial=FALSE,rug=FALSE,alpha=0.16)
visreg(gam2$gam,xvar="Age",by="Test",overlay=TRUE,ylab="Response Time",partial=FALSE,rug=FALSE,alpha=0.16)
dev.off()


# graph practice effects

x <- read.csv("CNB_Longitudinal_Core_11February2020.csv")

x <- x[which(x$Time < 5),]

gam1 <- gamm(ACC_ar~Test+s(Time,by=Test),random=list(bblid=~1),data=x)
gam2 <- gamm(RT_ar~Test+s(Time,by=Test),random=list(bblid=~1),data=x)

gam1$gam$data <- x
gam2$gam$data <- x

pdf("CNB_repeated-measures_Practice_Effects_by_test.pdf",width=12,height=6)
visreg(gam1$gam,xvar="Time",by="Test",overlay=TRUE,ylab="Accuracy",partial=FALSE,rug=FALSE,alpha=0.16)
visreg(gam2$gam,xvar="Time",by="Test",overlay=TRUE,ylab="Response Time",partial=FALSE,rug=FALSE,alpha=0.16)
dev.off()


# graph practice effects

x <- read.csv("CNB_Longitudinal_Core_11February2020.csv")

gam1 <- gamm(raw_original_acc~timepoint+int+timepoint:int+s(Age),random=list(bblid=~1),data=x[which(x$Test == "PCET"),])
gam2 <- gamm(raw_original_acc~timepoint+int+timepoint:int+s(Age),random=list(bblid=~1),data=x[which(x$Test == "CPT"),])
gam3 <- gamm(raw_original_acc~timepoint+int+timepoint:int+s(Age),random=list(bblid=~1),data=x[which(x$Test == "NBACK"),])
gam4 <- gamm(raw_original_acc~timepoint+int+timepoint:int+s(Age),random=list(bblid=~1),data=x[which(x$Test == "CPW"),])
gam5 <- gamm(raw_original_acc~timepoint+int+timepoint:int+s(Age),random=list(bblid=~1),data=x[which(x$Test == "CPF"),])
gam6 <- gamm(raw_original_acc~timepoint+int+timepoint:int+s(Age),random=list(bblid=~1),data=x[which(x$Test == "VOLT"),])
gam7 <- gamm(raw_original_acc~timepoint+int+timepoint:int+s(Age),random=list(bblid=~1),data=x[which(x$Test == "PVRT"),])
gam8 <- gamm(raw_original_acc~timepoint+int+timepoint:int+s(Age),random=list(bblid=~1),data=x[which(x$Test == "PMAT"),])
gam9 <- gamm(raw_original_acc~timepoint+int+timepoint:int+s(Age),random=list(bblid=~1),data=x[which(x$Test == "PLOT"),])
gam10 <- gamm(raw_original_acc~timepoint+int+timepoint:int+s(Age),random=list(bblid=~1),data=x[which(x$Test == "ER40"),])
gam11 <- gamm(raw_original_acc~timepoint+int+timepoint:int+s(Age),random=list(bblid=~1),data=x[which(x$Test == "MEDF"),])
gam12 <- gamm(raw_original_acc~timepoint+int+timepoint:int+s(Age),random=list(bblid=~1),data=x[which(x$Test == "ADT"),])
gam13 <- gamm(raw_original_acc~timepoint+int+timepoint:int+s(Age),random=list(bblid=~1),data=x[which(x$Test == "TAP"),])

gam1$gam$data <- x[which(x$Test == "PCET"),]
gam2$gam$data <- x[which(x$Test == "CPT"),]
gam3$gam$data <- x[which(x$Test == "NBACK"),]
gam4$gam$data <- x[which(x$Test == "CPW"),]
gam5$gam$data <- x[which(x$Test == "CPF"),]
gam6$gam$data <- x[which(x$Test == "VOLT"),]
gam7$gam$data <- x[which(x$Test == "PVRT"),]
gam8$gam$data <- x[which(x$Test == "PMAT"),]
gam9$gam$data <- x[which(x$Test == "PLOT"),]
gam10$gam$data <- x[which(x$Test == "ER40"),]
gam11$gam$data <- x[which(x$Test == "MEDF"),]
gam12$gam$data <- x[which(x$Test == "ADT"),]
gam13$gam$data <- x[which(x$Test == "TAP"),]

pdf("CNB_repeated-measures_Practice_Effects_by_test.pdf",width=12,height=6)
visreg(gam1$gam,xvar="Time",by="Test",overlay=TRUE,ylab="Accuracy",partial=FALSE,rug=FALSE,alpha=0.16)
visreg(gam2$gam,xvar="Time",by="Test",overlay=TRUE,ylab="Response Time",partial=FALSE,rug=FALSE,alpha=0.16)
dev.off()

