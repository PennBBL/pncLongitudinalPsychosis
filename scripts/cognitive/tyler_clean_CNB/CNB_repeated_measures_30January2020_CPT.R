
library(psych)
library(Amelia)
library(ggplot2)
library(gridExtra)
library(CorrMixed)
library(gtools)
library(mgcv)
library(visreg)
library(lubridate)
library(irr)

x <- read.csv("CNB Longitudinal/cpt_20191113.csv")

x[,4] <- as.numeric(ymd(x[,4])) # date of test
colnames(x)[4] <- "Time"
colnames(x)[3] <- "bblid"
colnames(x)[11] <- "Sex"

ids <- unique(x$bblid)
for (i in 1:length(ids)) {
temp <- x[which(x$bblid == ids[i]),]
if (is.element("N",temp[,6])) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])}
try(if (min(temp$cpt_ptp,na.rm=TRUE) < 0.25) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])})
try(if (max(temp$cpt_pfp,na.rm=TRUE) > 0.5) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])})
x[which(x$bblid == ids[i]),] <- temp}  

x <- x[which(x$ntimepoints > 1),]   # best is > 1
x <- x[which(x$timepoint < 6),]     # best is < 6
x <- x[which(is.na(x$cnbAgemonths) == FALSE),]

Age <- scale(x$cnbAgemonths)
Age_Squared <- Age^2
Age_Cubed <- Age^3
TP <- scale(x$timepoint)
TP_Squared <- scale(x$timepoint)^2
TP_Cubed <- scale(x$timepoint)^3

set.seed(2)
temp <- amelia(x[,16:19], m=1)$imputations[[1]]
x[,16:19] <- temp

dprime <- qnorm(x$cpt_ptp-0.001) - qnorm(x$cpt_pfp+0.001)

x[,16:19] <- scale(x[,16:19])

PTP_r <- scale(winsor(lm(cpt_ptp~Age+Age_Squared+Age_Cubed,data=x)$residuals,trim=0.005))
dmod <- lm(dprime~Age+Age_Squared+Age_Cubed,data=x,na.action=na.exclude)
dprime_r <- scale(winsor(residuals(dmod,na.action=na.exclude),trim=0.005))
TPRT_r <- scale(winsor(lm(cpt_tprt~Age+Age_Squared+Age_Cubed,data=x)$residuals,trim=0.005))
FPRT_r <- scale(winsor(lm(cpt_fprt~Age+Age_Squared+Age_Cubed,data=x)$residuals,trim=0.005))
RT_r <- scale(winsor(lm(cpt_tprt~Age+Age_Squared+Age_Cubed,data=x)$residuals,trim=0.005))

x <- data.frame(x,PTP_r,dprime_r,TPRT_r,FPRT_r,RT_r)

# arranges times to all start at 0
ids <- unique(x$bblid)
for (i in 1:length(ids)) {
temp <- x[which(x$bblid == ids[i]),]
temp$Time <- temp$Time - min(temp$Time,na.rm=TRUE)
x[which(x$bblid == ids[i]),] <- temp}

int <- matrix(NA,dim(x)[1],1)
x <- data.frame(x,int)

for (i in 1:length(ids)) {
temp <- x[which(x$bblid == ids[i]),]
for (j in 1:dim(temp)[1]) {
if (temp[j,20] == 9) {try(temp[j,27] <- temp[j,4] - temp[which(temp$timepoint == 8),4])}
if (temp[j,20] == 8) {try(temp[j,27] <- temp[j,4] - temp[which(temp$timepoint == 7),4])}
if (temp[j,20] == 7) {try(temp[j,27] <- temp[j,4] - temp[which(temp$timepoint == 6),4])}
if (temp[j,20] == 6) {try(temp[j,27] <- temp[j,4] - temp[which(temp$timepoint == 5),4])}
if (temp[j,20] == 5) {try(temp[j,27] <- temp[j,4] - temp[which(temp$timepoint == 4),4])}
if (temp[j,20] == 4) {try(temp[j,27] <- temp[j,4] - temp[which(temp$timepoint == 3),4])}
if (temp[j,20] == 3) {try(temp[j,27] <- temp[j,4] - temp[which(temp$timepoint == 2),4])}
if (temp[j,20] == 2) {try(temp[j,27] <- temp[j,4])}
if (temp[j,20] == 1) {try(temp[j,27] <- NA)}
x[which(x$bblid == ids[i]),] <- temp}}

x$cnbAgemonths <- x$cnbAgemonths/12
x$Time <- x$Time/365.25
Time_Squared <- scale(x$Time)^2
Age <- x$cnbAgemonths
int_ord <- as.factor(quantcut(x$int,3))

x <- data.frame(x,Time_Squared)

# raw basic models
gam1 <- gamm(cpt_ptp~Sex+s(Age,by=Sex),random=list(bblid=~1),data=x)
gam2 <- gamm(dprime~Sex+s(Age,by=Sex),random=list(bblid=~1),data=x)
gam3 <- gamm(cpt_tprt~Sex+s(Age,by=Sex),random=list(bblid=~1),data=x)

# models after age is regressed out cross-sectionally - give idea of practice effect
gam4 <- gamm(PTP_r~Sex+s(Age,by=Sex),random=list(bblid=~1),data=x)
gam5 <- gamm(dprime_r~Sex+s(Age,by=Sex),random=list(bblid=~1),data=x)
gam6 <- gamm(RT_r~Sex+s(Age,by=Sex),random=list(bblid=~1),data=x)

# how does inter-test interval affect score?
gam7 <- gamm(PTP_r~Sex+s(int,by=Sex),random=list(bblid=~1),data=x)
gam8 <- gamm(dprime_r~Sex+s(int,by=Sex),random=list(bblid=~1),data=x)
gam9 <- gamm(RT_r~Sex+s(int,by=Sex),random=list(bblid=~1),data=x)

# how does interval affect the practice effect?
gam10 <- gamm(PTP_r~int_ord+s(Time,by=int_ord),random=list(bblid=~1),data=x)
gam11 <- gamm(dprime_r~int_ord+s(Time,by=int_ord),random=list(bblid=~1),data=x)
gam12 <- gamm(RT_r~int_ord+s(Time,by=int_ord),random=list(bblid=~1),data=x)

# models after age is regressed out cross-sectionally - give idea of practice effect for TIME - this ignores practice/age interactions
gam13 <- gamm(PTP_r~Sex+s(Time,by=Sex),random=list(bblid=~1),data=x)
gam14 <- gamm(dprime_r~Sex+s(Time,by=Sex),random=list(bblid=~1),data=x)
gam15 <- gamm(RT_r~Sex+s(Time,by=Sex),random=list(bblid=~1),data=x)

# MORE MODELS ARE MADE WITHIN THE pdf() CREATION SECTION

gam1$gam$data <- x
gam2$gam$data <- x
gam3$gam$data <- x
gam4$gam$data <- x
gam5$gam$data <- x
gam6$gam$data <- x
gam7$gam$data <- x
gam8$gam$data <- x
gam9$gam$data <- x
gam10$gam$data <- x
gam11$gam$data <- x
gam12$gam$data <- x
gam13$gam$data <- x
gam14$gam$data <- x
gam15$gam$data <- x

# old code for visually examining models
#summary(lme(cpt_ptp~Age_Squared+(Sex+cnbAgemonths)^2,data=x,random = ~ 1 | bblid))
#summary(lme(dprime~Age_Squared+(Sex+cnbAgemonths)^2,data=x,random = ~ 1 | bblid))
#summary(lme(cpt_tprt~Age_Squared+(Sex+cnbAgemonths)^2,data=x,random = ~ 1 | bblid))

#summary(lme(PTP_r~Age_Squared+(Sex+Age)^2,data=x,random = ~ 1 | bblid))
#summary(lme(dprime_r~Age_Squared+(Sex+Age)^2,data=x,random = ~ 1 | bblid))
#summary(lme(RT_r~Age_Squared+(Sex+Age)^2,data=x,random = ~ 1 | bblid))

pdf("Repeated-Measures_Visuals_CPT.pdf",width=10,height=6)
visreg(gam1$gam,xvar="Age",by="Sex",overlay=TRUE,ylab="True Positives",partial=FALSE,rug=FALSE) #main interest
visreg(gam2$gam,xvar="Age",by="Sex",overlay=TRUE,ylab="dprime",partial=FALSE,rug=FALSE) #main interest
visreg(gam3$gam,xvar="Age",by="Sex",overlay=TRUE,ylab="Response Time",partial=FALSE,rug=FALSE)  #main interest
visreg(gam4$gam,xvar="Age",by="Sex",overlay=TRUE,ylab="True Positives (age-regressed)",partial=FALSE,rug=FALSE)  # rel'n with age after age-regress
visreg(gam5$gam,xvar="Age",by="Sex",overlay=TRUE,ylab="dprime (age-regressed)",partial=FALSE,rug=FALSE)  # rel'n with age after age-regress
visreg(gam6$gam,xvar="Age",by="Sex",overlay=TRUE,ylab="Response Time (age-regressed)",partial=FALSE,rug=FALSE)   # rel'n with age after age-regress
visreg(gam13$gam,xvar="Time",by="Sex",overlay=TRUE,ylab="True Positives (age-regressed)",partial=FALSE,rug=FALSE,strip.names=TRUE) # rel'n with time after age-regress
visreg(gam14$gam,xvar="Time",by="Sex",overlay=TRUE,ylab="dprime (age-regressed)",partial=FALSE,rug=FALSE,strip.names=TRUE) # rel'n with time after age-regress
visreg(gam15$gam,xvar="Time",by="Sex",overlay=TRUE,ylab="Response Time (age-regressed)",partial=FALSE,rug=FALSE,strip.names=TRUE) # rel'n with time after age-regress
visreg(gam7$gam,xvar="int",by="Sex",overlay=TRUE,ylab="True Positives (age-regressed)",xlab="Time Since Previous Administration",partial=FALSE,rug=FALSE) # interval
visreg(gam8$gam,xvar="int",by="Sex",overlay=TRUE,ylab="dprime (age-regressed)",xlab="Time Since Previous Administration",partial=FALSE,rug=FALSE) # interval
visreg(gam9$gam,xvar="int",by="Sex",overlay=TRUE,ylab="Response Time (age-regressed)",xlab="Time Since Previous Administration",partial=FALSE,rug=FALSE) # interval
visreg(gam10$gam,xvar="Time",by="int_ord",overlay=TRUE,ylab="True Positives (age-regressed)",partial=FALSE,rug=FALSE,strip.names=TRUE) # how does practice vary by interval?
visreg(gam11$gam,xvar="Time",by="int_ord",overlay=TRUE,ylab="dprime (age-regressed)",partial=FALSE,rug=FALSE,strip.names=TRUE) # how does practice vary by interval?
visreg(gam12$gam,xvar="Time",by="int_ord",overlay=TRUE,ylab="Response Time (age-regressed)",partial=FALSE,rug=FALSE,strip.names=TRUE) # how does practice vary by interval?

# new models built within visreg() - these are only to see the linear effects
visreg(lme(cpt_ptp~(Sex+cnbAgemonths)^2,data=x,random = ~ 1 | bblid),xvar="cnbAgemonths",by="Sex",overlay=TRUE,xlab="Age (Years)",ylab="True Positives")  # linear
visreg(lme(dprime~(Sex+cnbAgemonths)^2,data=x,random = ~ 1 | bblid),xvar="cnbAgemonths",by="Sex",overlay=TRUE,xlab="Age (Years)",ylab="dprime")  # linear
visreg(lme(cpt_tprt~(Sex+cnbAgemonths)^2,data=x,random = ~ 1 | bblid),xvar="cnbAgemonths",by="Sex",overlay=TRUE,xlab="Age (Years)",ylab="Response Time")  # linear
visreg(lme(PTP_r~(Sex+cnbAgemonths)^2,data=x,random = ~ 1 | bblid),xvar="cnbAgemonths",by="Sex",overlay=TRUE,xlab="Age (Years)",ylab="True Positives (Age-Regressed)") # linear
visreg(lme(dprime_r~(Sex+cnbAgemonths)^2,data=x,random = ~ 1 | bblid),xvar="cnbAgemonths",by="Sex",overlay=TRUE,xlab="Age (Years)",ylab="dprime (Age-Regressed)") # linear
visreg(lme(RT_r~(Sex+cnbAgemonths)^2,data=x,random = ~ 1 | bblid),xvar="cnbAgemonths",by="Sex",overlay=TRUE,xlab="Age (Years)",ylab="Response Time (Age-Regressed)") # linear
visreg(lme(PTP_r~(Sex+int)^2,data=x,random = ~ 1 | bblid,na.action=na.exclude),xvar="int",by="Sex",overlay=TRUE,xlab="Time Since Last Administration",ylab="True Positives (Age-Regressed)")  # interval
visreg(lme(dprime_r~(Sex+int)^2,data=x,random = ~ 1 | bblid,na.action=na.exclude),xvar="int",by="Sex",overlay=TRUE,xlab="Time Since Last Administration",ylab="dprime (Age-Regressed)")  # interval
visreg(lme(RT_r~(Sex+int)^2,data=x,random = ~ 1 | bblid,na.action=na.exclude),xvar="int",by="Sex",overlay=TRUE,xlab="Time Since Last Administration",ylab="Response Time (Age-Regressed)")  # interval

# below show actual mixed model results
plot.new()
title("dependent variable = cpt_ptp (raw true positives)")
grid.table(summary(lme(cpt_ptp~(Sex+cnbAgemonths)^2,data=x,random = ~ 1 | bblid))$tTable)
plot.new()
title("dependent variable = dprime")
grid.table(summary(lme(dprime~(Sex+cnbAgemonths)^2,data=x,random = ~ 1 | bblid))$tTable)
plot.new()
title("dependent variable = cpt_tprt (mean raw response time)")
grid.table(summary(lme(cpt_tprt~(Sex+cnbAgemonths)^2,data=x,random = ~ 1 | bblid))$tTable)
plot.new()
title("dependent variable = true positives (age-regressed)")
grid.table(summary(lme(PTP_r~(Sex+cnbAgemonths)^2,data=x,random = ~ 1 | bblid))$tTable)
plot.new()
title("dependent variable = dprime (age-regressed)")
grid.table(summary(lme(dprime_r~(Sex+cnbAgemonths)^2,data=x,random = ~ 1 | bblid))$tTable)
plot.new()
title("dependent variable = response time (age-regressed)")
grid.table(summary(lme(RT_r~(Sex+cnbAgemonths)^2,data=x,random = ~ 1 | bblid))$tTable)

# below we need people at time point 1 to have an interval score (set at highest in data set: 3128)
x$int[is.na(x$int) == TRUE] <- 3128
plot.new()
title("dependent variable = cpt_ptp (raw true positives)")
grid.table(summary(lme(cpt_ptp~(int+Sex+cnbAgemonths)^2,data=x,random = ~ 1 | bblid,na.action=na.exclude))$tTable) # try including interval
plot.new()
title("dependent variable = dprime")
grid.table(summary(lme(dprime~(int+Sex+cnbAgemonths)^2,data=x,random = ~ 1 | bblid,na.action=na.exclude))$tTable) # try including interval
plot.new()
title("dependent variable = cpt_tprt (mean raw response time)")
grid.table(summary(lme(cpt_tprt~(int+Sex+cnbAgemonths)^2,data=x,random = ~ 1 | bblid,na.action=na.exclude))$tTable) # try including interval

# these are the main models for practice effects - how do the number of previous admins and inter-admin interval affect age-regressed scores?
plot.new()
title("dependent variable = true positives (age-regressed)")
grid.table(summary(lme(PTP_r~(int+timepoint)^2,data=x,random = ~ 1 | bblid,na.action=na.exclude))$tTable) # big main model: age plus interval plus prev admins
plot.new()
title("dependent variable = dprime (age-regressed)")
grid.table(summary(lme(dprime_r~(int+timepoint)^2,data=x,random = ~ 1 | bblid,na.action=na.exclude))$tTable) # big main model: age plus interval plus prev admins
plot.new()
title("dependent variable = response time (age-regressed)")
grid.table(summary(lme(RT_r~(int+timepoint)^2,data=x,random = ~ 1 | bblid,na.action=na.exclude))$tTable) # big main model: age plus interval plus prev admins

# building these models here because if we built them up top with the other models, no one would have an interval score for time point 1
mod13 <- lme(PTP_r~(int+timepoint)^2,data=x,random = ~ 1 | bblid,na.action=na.exclude)
mod14 <- lme(dprime_r~(int+timepoint)^2,data=x,random = ~ 1 | bblid,na.action=na.exclude)
mod15 <- lme(RT_r~(int+timepoint)^2,data=x,random = ~ 1 | bblid,na.action=na.exclude)

# mod16 through mod18 are without the time point 1 interval
x$timepoint <- as.factor(x$timepoint)
mod16 <- lme(PTP_r~(int+timepoint)^2,data=x[which(x$timepoint != "1"),],random = ~ 1 | bblid,na.action=na.exclude)
mod17 <- lme(dprime_r~(int+timepoint)^2,data=x[which(x$timepoint != "1"),],random = ~ 1 | bblid,na.action=na.exclude)
mod18 <- lme(RT_r~(int+timepoint)^2,data=x[which(x$timepoint != "1"),],random = ~ 1 | bblid,na.action=na.exclude)

# plot the above six models in various ways.
visreg(mod13,xvar="int",by="timepoint",overlay=TRUE,xlab="Time Since Last Administration",ylab="True Positives (Age-Regressed)")
visreg(mod14,xvar="int",by="timepoint",overlay=TRUE,xlab="Time Since Last Administration",ylab="dprime (Age-Regressed)")
visreg(mod15,xvar="int",by="timepoint",overlay=TRUE,xlab="Time Since Last Administration",ylab="Response Time (Age-Regressed)")
visreg(mod13,xvar="timepoint",by="int",overlay=TRUE,xlab="Timepoint (administration number)",ylab="True Positives (Age-Regressed)")
visreg(mod14,xvar="timepoint",by="int",overlay=TRUE,xlab="Timepoint (administration number)",ylab="dprime (Age-Regressed)")
visreg(mod15,xvar="timepoint",by="int",overlay=TRUE,xlab="Timepoint (administration number)",ylab="Response Time (Age-Regressed)")
visreg(mod16,xvar="int",by="timepoint",overlay=TRUE,xlab="Time Since Last Administration",ylab="True Positives (Age-Regressed)")
visreg(mod17,xvar="int",by="timepoint",overlay=TRUE,xlab="Time Since Last Administration",ylab="dprime (Age-Regressed)")
visreg(mod18,xvar="int",by="timepoint",overlay=TRUE,xlab="Time Since Last Administration",ylab="Response Time (Age-Regressed)")
visreg(mod16,xvar="timepoint",by="int",overlay=TRUE,xlab="Timepoint (administration number)",ylab="True Positives (Age-Regressed)")
visreg(mod17,xvar="timepoint",by="int",overlay=TRUE,xlab="Timepoint (administration number)",ylab="dprime (Age-Regressed)")
visreg(mod18,xvar="timepoint",by="int",overlay=TRUE,xlab="Timepoint (administration number)",ylab="Response Time (Age-Regressed)")
dev.off()



# make wide and run test-retest stats

x <- read.csv("CNB Longitudinal/cpt_20191113.csv")

x[,4] <- as.numeric(ymd(x[,4])) # date of test
colnames(x)[4] <- "Time"
colnames(x)[3] <- "bblid"
colnames(x)[11] <- "Sex"

ids <- unique(x$bblid)
for (i in 1:length(ids)) {
temp <- x[which(x$bblid == ids[i]),]
if (is.element("N",temp[,6])) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])}
try(if (min(temp$cpt_ptp,na.rm=TRUE) < 0.25) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])})
try(if (max(temp$cpt_pfp,na.rm=TRUE) > 0.5) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])})
x[which(x$bblid == ids[i]),] <- temp}  

x <- x[which(x$ntimepoints > 1),]   # best is > 1
x <- x[which(x$timepoint < 8),]     # best is < 6
x <- x[which(is.na(x$cnbAgemonths) == FALSE),]

Age <- scale(x$cnbAgemonths)
Age_Squared <- Age^2
Age_Cubed <- Age^3
TP <- scale(x$timepoint)
TP_Squared <- scale(x$timepoint)^2
TP_Cubed <- scale(x$timepoint)^3

set.seed(2)
temp <- amelia(x[,16:19], m=1)$imputations[[1]]
x[,16:19] <- temp

dprime <- qnorm(x$cpt_ptp-0.001) - qnorm(x$cpt_pfp+0.001)

x[,16:19] <- scale(x[,16:19])

# arranges times to all start at 0
ids <- unique(x$bblid)
for (i in 1:length(ids)) {
temp <- x[which(x$bblid == ids[i]),]
temp$Time <- temp$Time - min(temp$Time,na.rm=TRUE)
x[which(x$bblid == ids[i]),] <- temp}

int <- matrix(NA,dim(x)[1],1)
x <- data.frame(x,int)

for (i in 1:length(ids)) {
temp <- x[which(x$bblid == ids[i]),]
for (j in 1:dim(temp)[1]) {
if (temp[j,20] == 9) {try(temp[j,22] <- temp[j,4] - temp[which(temp$timepoint == 8),4])}
if (temp[j,20] == 8) {try(temp[j,22] <- temp[j,4] - temp[which(temp$timepoint == 7),4])}
if (temp[j,20] == 7) {try(temp[j,22] <- temp[j,4] - temp[which(temp$timepoint == 6),4])}
if (temp[j,20] == 6) {try(temp[j,22] <- temp[j,4] - temp[which(temp$timepoint == 5),4])}
if (temp[j,20] == 5) {try(temp[j,22] <- temp[j,4] - temp[which(temp$timepoint == 4),4])}
if (temp[j,20] == 4) {try(temp[j,22] <- temp[j,4] - temp[which(temp$timepoint == 3),4])}
if (temp[j,20] == 3) {try(temp[j,22] <- temp[j,4] - temp[which(temp$timepoint == 2),4])}
if (temp[j,20] == 2) {try(temp[j,22] <- temp[j,4])}
if (temp[j,20] == 1) {try(temp[j,22] <- 3128)}
x[which(x$bblid == ids[i]),] <- temp}}

int_sq <- scale(x[,22])^2
int_cub <- scale(x[,22])^3

#regressing out age
PTP_ar <- scale(winsor(lm(cpt_ptp~Age+Age_Squared+Age_Cubed,data=x)$residuals,trim=0.005))
dmod <- lm(dprime~Age+Age_Squared+Age_Cubed,data=x,na.action=na.exclude)
dprime_ar <- scale(winsor(residuals(dmod,na.action=na.exclude),trim=0.005))
TPRT_ar <- scale(winsor(lm(cpt_tprt~Age+Age_Squared+Age_Cubed,data=x)$residuals,trim=0.005))
FPRT_ar <- scale(winsor(lm(cpt_fprt~Age+Age_Squared+Age_Cubed,data=x)$residuals,trim=0.005))
RT_ar <- scale(winsor(lm(cpt_tprt~Age+Age_Squared+Age_Cubed,data=x)$residuals,trim=0.005))

#regressing out age and practice
PTP_apr <- scale(winsor(lm(cpt_ptp~Age+Age_Squared+Age_Cubed+TP+TP_Squared+TP_Cubed+int+int_sq+int_cub,data=x)$residuals,trim=0.005))
dmod <- lm(dprime~Age+Age_Squared+Age_Cubed+TP+TP_Squared+TP_Cubed+int+int_sq+int_cub,data=x,na.action=na.exclude)
dprime_apr <- scale(winsor(residuals(dmod,na.action=na.exclude),trim=0.005))
TPRT_apr <- scale(winsor(lm(cpt_tprt~Age+Age_Squared+Age_Cubed+TP+TP_Squared+TP_Cubed+int+int_sq+int_cub,data=x)$residuals,trim=0.005))
FPRT_apr <- scale(winsor(lm(cpt_fprt~Age+Age_Squared+Age_Cubed+TP+TP_Squared+TP_Cubed+int+int_sq+int_cub,data=x)$residuals,trim=0.005))
RT_apr <- scale(winsor(lm(cpt_tprt~Age+Age_Squared+Age_Cubed+TP+TP_Squared+TP_Cubed+int+int_sq+int_cub,data=x)$residuals,trim=0.005))

x <- data.frame(x,PTP_ar,dprime_ar,TPRT_ar,FPRT_ar,RT_ar,PTP_apr,dprime_apr,TPRT_apr,FPRT_apr,RT_apr)

x1 <- x[which(x$timepoint == 1),]
x2 <- x[which(x$timepoint == 2),]
x3 <- x[which(x$timepoint == 3),]
x4 <- x[which(x$timepoint == 4),]
x5 <- x[which(x$timepoint == 5),]
x6 <- x[which(x$timepoint == 6),]
x7 <- x[which(x$timepoint == 7),]

x <- merge(x1,x2,by="bblid",all=TRUE)
x <- merge(x,x3,by="bblid",all=TRUE)
x <- merge(x,x4,by="bblid",all=TRUE)
x <- merge(x,x5,by="bblid",all=TRUE)
x <- merge(x,x6,by="bblid",all=TRUE)
x <- merge(x,x7,by="bblid",all=TRUE)

#create columns names because there are duplicate columns in current frame
colnames(x) <- c(
"bblid",
"datasetid_TIME_1",
"siteid_TIME_1",
"Time_TIME_1",
"battery_TIME_1",
"valid_code_TIME_1",
"age_TIME_1",
"education_TIME_1",
"feducation_TIME_1",
"meducation_TIME_1",
"Sex_TIME_1",
"handedness_TIME_1",
"cnbAgemonths_TIME_1",
"cpt_genus_TIME_1",
"cpt_valid_TIME_1",
"cpt_ptp_TIME_1",
"cpt_tprt_TIME_1",
"cpt_pfp_TIME_1",
"cpt_fprt_TIME_1",
"timepoint_TIME_1",
"ntimepoints_TIME_1",
"int_TIME_1",
"PTP_ar_TIME_1",
"dprime_ar_TIME_1",
"TPRT_ar_TIME_1",
"FPRT_ar_TIME_1",
"RT_ar_TIME_1",
"PTP_apr_TIME_1",
"dprime_apr_TIME_1",
"TPRT_apr_TIME_1",
"FPRT_apr_TIME_1",
"RT_apr_TIME_1",
"datasetid_TIME_2",
"siteid_TIME_2",
"Time_TIME_2",
"battery_TIME_2",
"valid_code_TIME_2",
"age_TIME_2",
"education_TIME_2",
"feducation_TIME_2",
"meducation_TIME_2",
"Sex_TIME_2",
"handedness_TIME_2",
"cnbAgemonths_TIME_2",
"cpt_genus_TIME_2",
"cpt_valid_TIME_2",
"cpt_ptp_TIME_2",
"cpt_tprt_TIME_2",
"cpt_pfp_TIME_2",
"cpt_fprt_TIME_2",
"timepoint_TIME_2",
"ntimepoints_TIME_2",
"int_TIME_2",
"PTP_ar_TIME_2",
"dprime_ar_TIME_2",
"TPRT_ar_TIME_2",
"FPRT_ar_TIME_2",
"RT_ar_TIME_2",
"PTP_apr_TIME_2",
"dprime_apr_TIME_2",
"TPRT_apr_TIME_2",
"FPRT_apr_TIME_2",
"RT_apr_TIME_2",
"datasetid_TIME_3",
"siteid_TIME_3",
"Time_TIME_3",
"battery_TIME_3",
"valid_code_TIME_3",
"age_TIME_3",
"education_TIME_3",
"feducation_TIME_3",
"meducation_TIME_3",
"Sex_TIME_3",
"handedness_TIME_3",
"cnbAgemonths_TIME_3",
"cpt_genus_TIME_3",
"cpt_valid_TIME_3",
"cpt_ptp_TIME_3",
"cpt_tprt_TIME_3",
"cpt_pfp_TIME_3",
"cpt_fprt_TIME_3",
"timepoint_TIME_3",
"ntimepoints_TIME_3",
"int_TIME_3",
"PTP_ar_TIME_3",
"dprime_ar_TIME_3",
"TPRT_ar_TIME_3",
"FPRT_ar_TIME_3",
"RT_ar_TIME_3",
"PTP_apr_TIME_3",
"dprime_apr_TIME_3",
"TPRT_apr_TIME_3",
"FPRT_apr_TIME_3",
"RT_apr_TIME_3",
"datasetid_TIME_4",
"siteid_TIME_4",
"Time_TIME_4",
"battery_TIME_4",
"valid_code_TIME_4",
"age_TIME_4",
"education_TIME_4",
"feducation_TIME_4",
"meducation_TIME_4",
"Sex_TIME_4",
"handedness_TIME_4",
"cnbAgemonths_TIME_4",
"cpt_genus_TIME_4",
"cpt_valid_TIME_4",
"cpt_ptp_TIME_4",
"cpt_tprt_TIME_4",
"cpt_pfp_TIME_4",
"cpt_fprt_TIME_4",
"timepoint_TIME_4",
"ntimepoints_TIME_4",
"int_TIME_4",
"PTP_ar_TIME_4",
"dprime_ar_TIME_4",
"TPRT_ar_TIME_4",
"FPRT_ar_TIME_4",
"RT_ar_TIME_4",
"PTP_apr_TIME_4",
"dprime_apr_TIME_4",
"TPRT_apr_TIME_4",
"FPRT_apr_TIME_4",
"RT_apr_TIME_4",
"datasetid_TIME_5",
"siteid_TIME_5",
"Time_TIME_5",
"battery_TIME_5",
"valid_code_TIME_5",
"age_TIME_5",
"education_TIME_5",
"feducation_TIME_5",
"meducation_TIME_5",
"Sex_TIME_5",
"handedness_TIME_5",
"cnbAgemonths_TIME_5",
"cpt_genus_TIME_5",
"cpt_valid_TIME_5",
"cpt_ptp_TIME_5",
"cpt_tprt_TIME_5",
"cpt_pfp_TIME_5",
"cpt_fprt_TIME_5",
"timepoint_TIME_5",
"ntimepoints_TIME_5",
"int_TIME_5",
"PTP_ar_TIME_5",
"dprime_ar_TIME_5",
"TPRT_ar_TIME_5",
"FPRT_ar_TIME_5",
"RT_ar_TIME_5",
"PTP_apr_TIME_5",
"dprime_apr_TIME_5",
"TPRT_apr_TIME_5",
"FPRT_apr_TIME_5",
"RT_apr_TIME_5",
"datasetid_TIME_6",
"siteid_TIME_6",
"Time_TIME_6",
"battery_TIME_6",
"valid_code_TIME_6",
"age_TIME_6",
"education_TIME_6",
"feducation_TIME_6",
"meducation_TIME_6",
"Sex_TIME_6",
"handedness_TIME_6",
"cnbAgemonths_TIME_6",
"cpt_genus_TIME_6",
"cpt_valid_TIME_6",
"cpt_ptp_TIME_6",
"cpt_tprt_TIME_6",
"cpt_pfp_TIME_6",
"cpt_fprt_TIME_6",
"timepoint_TIME_6",
"ntimepoints_TIME_6",
"int_TIME_6",
"PTP_ar_TIME_6",
"dprime_ar_TIME_6",
"TPRT_ar_TIME_6",
"FPRT_ar_TIME_6",
"RT_ar_TIME_6",
"PTP_apr_TIME_6",
"dprime_apr_TIME_6",
"TPRT_apr_TIME_6",
"FPRT_apr_TIME_6",
"RT_apr_TIME_6",
"datasetid_TIME_7",
"siteid_TIME_7",
"Time_TIME_7",
"battery_TIME_7",
"valid_code_TIME_7",
"age_TIME_7",
"education_TIME_7",
"feducation_TIME_7",
"meducation_TIME_7",
"Sex_TIME_7",
"handedness_TIME_7",
"cnbAgemonths_TIME_7",
"cpt_genus_TIME_7",
"cpt_valid_TIME_7",
"cpt_ptp_TIME_7",
"cpt_tprt_TIME_7",
"cpt_pfp_TIME_7",
"cpt_fprt_TIME_7",
"timepoint_TIME_7",
"ntimepoints_TIME_7",
"int_TIME_7",
"PTP_ar_TIME_7",
"dprime_ar_TIME_7",
"TPRT_ar_TIME_7",
"FPRT_ar_TIME_7",
"RT_ar_TIME_7",
"PTP_apr_TIME_7",
"dprime_apr_TIME_7",
"TPRT_apr_TIME_7",
"FPRT_apr_TIME_7",
"RT_apr_TIME_7")

PTP <- data.frame(
x$cpt_ptp_TIME_1,
x$cpt_ptp_TIME_2,
x$cpt_ptp_TIME_3,
x$cpt_ptp_TIME_4,
x$cpt_ptp_TIME_5,
x$cpt_ptp_TIME_6,
x$cpt_ptp_TIME_7)

PFP <- data.frame(
x$cpt_pfp_TIME_1,
x$cpt_pfp_TIME_2,
x$cpt_pfp_TIME_3,
x$cpt_pfp_TIME_4,
x$cpt_pfp_TIME_5,
x$cpt_pfp_TIME_6,
x$cpt_pfp_TIME_7)

RT <- data.frame(
x$cpt_tprt_TIME_1,
x$cpt_tprt_TIME_2,
x$cpt_tprt_TIME_3,
x$cpt_tprt_TIME_4,
x$cpt_tprt_TIME_5,
x$cpt_tprt_TIME_6,
x$cpt_tprt_TIME_7)

PTP_ar <- data.frame(
x$PTP_ar_TIME_1,
x$PTP_ar_TIME_2,
x$PTP_ar_TIME_3,
x$PTP_ar_TIME_4,
x$PTP_ar_TIME_5,
x$PTP_ar_TIME_6,
x$PTP_ar_TIME_7)

dprime_ar <- data.frame(
x$dprime_ar_TIME_1,
x$dprime_ar_TIME_2,
x$dprime_ar_TIME_3,
x$dprime_ar_TIME_4,
x$dprime_ar_TIME_5,
x$dprime_ar_TIME_6,
x$dprime_ar_TIME_7)

RT_ar <- data.frame(
x$RT_ar_TIME_1,
x$RT_ar_TIME_2,
x$RT_ar_TIME_3,
x$RT_ar_TIME_4,
x$RT_ar_TIME_5,
x$RT_ar_TIME_6,
x$RT_ar_TIME_7)

PTP_apr <- data.frame(
x$PTP_apr_TIME_1,
x$PTP_apr_TIME_2,
x$PTP_apr_TIME_3,
x$PTP_apr_TIME_4,
x$PTP_apr_TIME_5,
x$PTP_apr_TIME_6,
x$PTP_apr_TIME_7)

dprime_apr <- data.frame(
x$dprime_apr_TIME_1,
x$dprime_apr_TIME_2,
x$dprime_apr_TIME_3,
x$dprime_apr_TIME_4,
x$dprime_apr_TIME_5,
x$dprime_apr_TIME_6,
x$dprime_apr_TIME_7)

RT_apr <- data.frame(
x$RT_apr_TIME_1,
x$RT_apr_TIME_2,
x$RT_apr_TIME_3,
x$RT_apr_TIME_4,
x$RT_apr_TIME_5,
x$RT_apr_TIME_6,
x$RT_apr_TIME_7)

res <- matrix(NA,6,5)

for (i in 2:6) {
res[1,(i-1)] <- icc(PTP_apr[,1:i],type="agreement",model="twoway")$value
res[2,(i-1)] <- icc(dprime_apr[,1:i],type="agreement",model="twoway")$value
res[3,(i-1)] <- icc(RT_apr[,1:i],type="agreement",model="twoway")$value
res[4,(i-1)] <- icc(PTP_ar[,1:i],type="agreement",model="twoway")$value
res[5,(i-1)] <- icc(dprime_ar[,1:i],type="agreement",model="twoway")$value
res[6,(i-1)] <- icc(RT_ar[,1:i],type="agreement",model="twoway")$value
} 

header <- c("2 Timepoints","3 Timepoints","4 Timepoints","5 Timepoints","6 Timepoints")
TestRetest_for_CPT <- c("","True Positives (age- & practice-regressed)","dprime (age- & practice-regressed)","Response Time (age- & practice-regressed)","True Positives (age-regressed)","dprime (age-regressed)","Response Time (age-regressed)")

pdf("TestRetest_CPT.pdf",width=15,height=15)
grid.table(cbind(TestRetest_for_CPT,rbind(header,round(res,3))))
pairs.panels(PTP_apr,lm=TRUE)
pairs.panels(dprime_apr,lm=TRUE)
pairs.panels(RT_apr,lm=TRUE)
dev.off()




