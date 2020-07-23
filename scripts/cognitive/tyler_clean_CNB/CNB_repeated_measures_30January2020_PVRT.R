
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

x <- read.csv("CNB Longitudinal/pvrt_20191113.csv")

x[,4] <- as.numeric(ymd(x[,4])) # date of test
colnames(x)[4] <- "Time"
colnames(x)[3] <- "bblid"
colnames(x)[11] <- "Sex"

x <- x[grepl("KSPVRT",x$pvrt_genus),]

ids <- unique(x$bblid)
for (i in 1:length(ids)) {
temp <- x[which(x$bblid == ids[i]),]
if (is.element("N",temp[,6])) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])}
try(if (min(temp$pvrt_pc,na.rm=TRUE) < 10) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])})
try(if (max(temp$pvrt_rtcr,na.rm=TRUE) > 20000) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])})
x[which(x$bblid == ids[i]),] <- temp}  

x <- x[which(x$ntimepoints > 1),]   # best is > 1
x <- x[which(x$timepoint < 5),]     # best is < 6
x <- x[which(is.na(x$cnbAgemonths) == FALSE),]

ids <- unique(x$bblid)
for (i in 1:length(ids)) {
temp <- x[which(x$bblid == ids[i]),]
temp[,19] <- seq(1:nrow(temp))
x[which(x$bblid == ids[i]),] <- temp}

Age <- scale(x$cnbAgemonths)
Age_Squared <- Age^2
Age_Cubed <- Age^3
TP <- scale(x$timepoint)
TP_Squared <- scale(x$timepoint)^2
TP_Cubed <- scale(x$timepoint)^3

#set.seed(2)
#temp <- amelia(x[,16:19], m=1)$imputations[[1]]
#x[,16:19] <- temp

x[,16:18] <- scale(x[,16:18])

ACC_r <- scale(winsor(lm(pvrt_pc~Age+Age_Squared+Age_Cubed,data=x)$residuals,trim=0.005))
RT_r <- scale(winsor(lm(pvrt_rtcr~Age+Age_Squared+Age_Cubed,data=x)$residuals,trim=0.005))

x <- data.frame(x,ACC_r,RT_r)

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
if (temp[j,19] == 4) {try(temp[j,23] <- temp[j,4] - temp[which(temp$timepoint == 3),4])}
if (temp[j,19] == 3) {try(temp[j,23] <- temp[j,4] - temp[which(temp$timepoint == 2),4])}
if (temp[j,19] == 2) {try(temp[j,23] <- temp[j,4])}
if (temp[j,19] == 1) {try(temp[j,23] <- NA)}
x[which(x$bblid == ids[i]),] <- temp}}

max_int <- max(x$int,na.rm=TRUE)

x$cnbAgemonths <- x$cnbAgemonths/12
x$Time <- x$Time/365.25
Time_Squared <- scale(x$Time)^2
Age <- x$cnbAgemonths
int_ord <- as.factor(quantcut(x$int,3))

x <- data.frame(x,Time_Squared)

# raw basic models
gam1 <- gamm(pvrt_pc~Sex+s(Age,by=Sex),random=list(bblid=~1),data=x)
gam2 <- gamm(pvrt_rtcr~Sex+s(Age,by=Sex),random=list(bblid=~1),data=x)

# models after age is regressed out cross-sectionally - give idea of practice effect
gam3 <- gamm(ACC_r~Sex+s(Age,by=Sex),random=list(bblid=~1),data=x)
gam4 <- gamm(RT_r~Sex+s(Age,by=Sex),random=list(bblid=~1),data=x)

# how does inter-test interval affect score?
gam5 <- gamm(ACC_r~Sex+s(int,by=Sex),random=list(bblid=~1),data=x)
gam6 <- gamm(RT_r~Sex+s(int,by=Sex),random=list(bblid=~1),data=x)

# how does interval affect the practice effect?
gam7 <- gamm(ACC_r~int_ord+s(Time,by=int_ord),random=list(bblid=~1),data=x)
gam8 <- gamm(RT_r~int_ord+s(Time,by=int_ord),random=list(bblid=~1),data=x)

# models after age is regressed out cross-sectionally - give idea of practice effect for TIME - this ignores practice/age interactions
gam9 <- gamm(ACC_r~Sex+s(Time,by=Sex),random=list(bblid=~1),data=x)
gam10 <- gamm(RT_r~Sex+s(Time,by=Sex),random=list(bblid=~1),data=x)

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

# old code for visually examining models
#summary(lme(pvrt_pc~Age_Squared+(Sex+cnbAgemonths)^2,data=x,random = ~ 1 | bblid))
#summary(lme(dprime~Age_Squared+(Sex+cnbAgemonths)^2,data=x,random = ~ 1 | bblid))
#summary(lme(pvrt_rtcr~Age_Squared+(Sex+cnbAgemonths)^2,data=x,random = ~ 1 | bblid))

#summary(lme(ACC_r~Age_Squared+(Sex+Age)^2,data=x,random = ~ 1 | bblid))
#summary(lme(dprime_r~Age_Squared+(Sex+Age)^2,data=x,random = ~ 1 | bblid))
#summary(lme(RT_r~Age_Squared+(Sex+Age)^2,data=x,random = ~ 1 | bblid))

pdf("Repeated-Measures_Visuals_PVRT.pdf",width=10,height=6)
visreg(gam1$gam,xvar="Age",by="Sex",overlay=TRUE,ylab="Accuracy",partial=FALSE,rug=FALSE) #main interest
visreg(gam2$gam,xvar="Age",by="Sex",overlay=TRUE,ylab="Response Time",partial=FALSE,rug=FALSE)  #main interest
visreg(gam3$gam,xvar="Age",by="Sex",overlay=TRUE,ylab="Accuracy (age-regressed)",partial=FALSE,rug=FALSE)  # rel'n with age after age-regress
visreg(gam4$gam,xvar="Age",by="Sex",overlay=TRUE,ylab="Response Time (age-regressed)",partial=FALSE,rug=FALSE)   # rel'n with age after age-regress
visreg(gam9$gam,xvar="Time",by="Sex",overlay=TRUE,ylab="Accuracy (age-regressed)",partial=FALSE,rug=FALSE,strip.names=TRUE) # rel'n with time after age-regress
visreg(gam10$gam,xvar="Time",by="Sex",overlay=TRUE,ylab="Response Time (age-regressed)",partial=FALSE,rug=FALSE,strip.names=TRUE) # rel'n with time after age-regress
visreg(gam5$gam,xvar="int",by="Sex",overlay=TRUE,ylab="Accuracy (age-regressed)",xlab="Time Since Previous Administration",partial=FALSE,rug=FALSE) # interval
visreg(gam6$gam,xvar="int",by="Sex",overlay=TRUE,ylab="Response Time (age-regressed)",xlab="Time Since Previous Administration",partial=FALSE,rug=FALSE) # interval
visreg(gam7$gam,xvar="Time",by="int_ord",overlay=TRUE,ylab="Accuracy (age-regressed)",partial=FALSE,rug=FALSE,strip.names=TRUE) # how does practice vary by interval?
visreg(gam8$gam,xvar="Time",by="int_ord",overlay=TRUE,ylab="Response Time (age-regressed)",partial=FALSE,rug=FALSE,strip.names=TRUE) # how does practice vary by interval?

# new models built within visreg() - these are only to see the linear effects
visreg(lme(pvrt_pc~(Sex+cnbAgemonths)^2,data=x,random = ~ 1 | bblid),xvar="cnbAgemonths",by="Sex",overlay=TRUE,xlab="Age (Years)",ylab="Accuracy")  # linear
visreg(lme(pvrt_rtcr~(Sex+cnbAgemonths)^2,data=x,random = ~ 1 | bblid),xvar="cnbAgemonths",by="Sex",overlay=TRUE,xlab="Age (Years)",ylab="Response Time")  # linear
visreg(lme(ACC_r~(Sex+cnbAgemonths)^2,data=x,random = ~ 1 | bblid),xvar="cnbAgemonths",by="Sex",overlay=TRUE,xlab="Age (Years)",ylab="Accuracy (Age-Regressed)") # linear
visreg(lme(RT_r~(Sex+cnbAgemonths)^2,data=x,random = ~ 1 | bblid),xvar="cnbAgemonths",by="Sex",overlay=TRUE,xlab="Age (Years)",ylab="Response Time (Age-Regressed)") # linear
visreg(lme(ACC_r~(Sex+int)^2,data=x,random = ~ 1 | bblid,na.action=na.exclude),xvar="int",by="Sex",overlay=TRUE,xlab="Time Since Last Administration",ylab="Accuracy (Age-Regressed)") # interval
visreg(lme(RT_r~(Sex+int)^2,data=x,random = ~ 1 | bblid,na.action=na.exclude),xvar="int",by="Sex",overlay=TRUE,xlab="Time Since Last Administration",ylab="Response Time (Age-Regressed)") # interval

# below show actual mixed model results
plot.new()
title("dependent variable = pvrt_pc (raw Accuracy)")
grid.table(summary(lme(pvrt_pc~(Sex+cnbAgemonths)^2,data=x,random = ~ 1 | bblid))$tTable)
plot.new()
title("dependent variable = pvrt_rtcr (mean raw response time)")
grid.table(summary(lme(pvrt_rtcr~(Sex+cnbAgemonths)^2,data=x,random = ~ 1 | bblid))$tTable)
plot.new()
title("dependent variable = Accuracy (age-regressed)")
grid.table(summary(lme(ACC_r~(Sex+cnbAgemonths)^2,data=x,random = ~ 1 | bblid))$tTable)
plot.new()
title("dependent variable = response time (age-regressed)")
grid.table(summary(lme(RT_r~(Sex+cnbAgemonths)^2,data=x,random = ~ 1 | bblid))$tTable)

# below we need people at time point 1 to have an interval score (set at highest in data set)
x$int[is.na(x$int) == TRUE] <- max_int
plot.new()
title("dependent variable = pvrt_pc (raw Accuracy)")
grid.table(summary(lme(pvrt_pc~(int+Sex+cnbAgemonths)^2,data=x,random = ~ 1 | bblid,na.action=na.exclude))$tTable) # try including interval
plot.new()
title("dependent variable = pvrt_rtcr (mean raw response time)")
grid.table(summary(lme(pvrt_rtcr~(int+Sex+cnbAgemonths)^2,data=x,random = ~ 1 | bblid,na.action=na.exclude))$tTable) # try including interval

# these are the main models for practice effects - how do the number of previous admins and inter-admin interval affect age-regressed scores?
plot.new()
title("dependent variable = Accuracy (age-regressed)")
grid.table(summary(lme(ACC_r~(int+timepoint)^2,data=x,random = ~ 1 | bblid,na.action=na.exclude))$tTable) # big main model: age plus interval plus prev admins
plot.new()
title("dependent variable = response time (age-regressed)")
grid.table(summary(lme(RT_r~(int+timepoint)^2,data=x,random = ~ 1 | bblid,na.action=na.exclude))$tTable) # big main model: age plus interval plus prev admins

# building these models here because if we built them up top with the other models, no one would have an interval score for time point 1
mod11 <- lme(ACC_r~(int+timepoint)^2,data=x,random = ~ 1 | bblid,na.action=na.exclude)
mod12 <- lme(RT_r~(int+timepoint)^2,data=x,random = ~ 1 | bblid,na.action=na.exclude)

# mod16 through mod18 are without the time point 1 interval
x$timepoint <- as.factor(x$timepoint)
mod13 <- lme(ACC_r~(int+timepoint)^2,data=x[which(x$timepoint != "1"),],random = ~ 1 | bblid,na.action=na.exclude)
mod14 <- lme(RT_r~(int+timepoint)^2,data=x[which(x$timepoint != "1"),],random = ~ 1 | bblid,na.action=na.exclude)

# plot the above six models in various ways.
visreg(mod11,xvar="int",by="timepoint",overlay=TRUE,xlab="Time Since Last Administration",ylab="Accuracy (Age-Regressed)",strip.names=TRUE)
visreg(mod12,xvar="int",by="timepoint",overlay=TRUE,xlab="Time Since Last Administration",ylab="Response Time (Age-Regressed)",strip.names=TRUE)
visreg(mod11,xvar="timepoint",by="int",overlay=TRUE,xlab="Timepoint (administration number)",ylab="Accuracy (Age-Regressed)",strip.names=TRUE)
visreg(mod12,xvar="timepoint",by="int",overlay=TRUE,xlab="Timepoint (administration number)",ylab="Response Time (Age-Regressed)",strip.names=TRUE)
visreg(mod13,xvar="int",by="timepoint",overlay=TRUE,xlab="Time Since Last Administration",ylab="Accuracy (Age-Regressed)",strip.names=TRUE)
visreg(mod14,xvar="int",by="timepoint",overlay=TRUE,xlab="Time Since Last Administration",ylab="Response Time (Age-Regressed)",strip.names=TRUE)
visreg(mod13,xvar="timepoint",by="int",overlay=TRUE,xlab="Timepoint (administration number)",ylab="Accuracy (Age-Regressed)",strip.names=TRUE)
visreg(mod14,xvar="timepoint",by="int",overlay=TRUE,xlab="Timepoint (administration number)",ylab="Response Time (Age-Regressed)",strip.names=TRUE)
dev.off()



# make wide and run test-retest stats

x <- read.csv("CNB Longitudinal/pvrt_20191113.csv")

x[,4] <- as.numeric(ymd(x[,4])) # date of test
colnames(x)[4] <- "Time"
colnames(x)[3] <- "bblid"
colnames(x)[11] <- "Sex"

ids <- unique(x$bblid)
for (i in 1:length(ids)) {
temp <- x[which(x$bblid == ids[i]),]
if (is.element("N",temp[,6])) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])}
try(if (min(temp$pvrt_pc,na.rm=TRUE) < 10) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])})
try(if (max(temp$pvrt_rtcr,na.rm=TRUE) > 20000) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])})
x[which(x$bblid == ids[i]),] <- temp}  

x <- x[which(x$ntimepoints > 1),]   # best is > 1
x <- x[which(x$timepoint < 6),]     # best is < 6
x <- x[which(is.na(x$cnbAgemonths) == FALSE),]

x <- x[grepl("KSPVRT",x$pvrt_genus),]

ids <- unique(x$bblid)
for (i in 1:length(ids)) {
temp <- x[which(x$bblid == ids[i]),]
temp[,19] <- seq(1:nrow(temp))
x[which(x$bblid == ids[i]),] <- temp}

Age <- scale(x$cnbAgemonths)
Age_Squared <- Age^2
Age_Cubed <- Age^3
TP <- scale(x$timepoint)
TP_Squared <- scale(x$timepoint)^2
TP_Cubed <- scale(x$timepoint)^3

#set.seed(2)
#temp <- amelia(x[,16:19], m=1)$imputations[[1]]
#x[,16:19] <- temp

x[,16:18] <- scale(x[,16:18])

#regressing out age
ACC_ar <- scale(winsor(lm(pvrt_pc~Age+Age_Squared+Age_Cubed,data=x)$residuals,trim=0.005))
RT_ar <- scale(winsor(lm(pvrt_rtcr~Age+Age_Squared+Age_Cubed,data=x)$residuals,trim=0.005))

x <- data.frame(x,ACC_ar,RT_ar)

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
if (temp[j,19] == 5) {try(temp[j,23] <- temp[j,4] - temp[which(temp$timepoint == 4),4])}
if (temp[j,19] == 4) {try(temp[j,23] <- temp[j,4] - temp[which(temp$timepoint == 3),4])}
if (temp[j,19] == 3) {try(temp[j,23] <- temp[j,4] - temp[which(temp$timepoint == 2),4])}
if (temp[j,19] == 2) {try(temp[j,23] <- temp[j,4])}
if (temp[j,19] == 1) {try(temp[j,23] <- max_int)}
x[which(x$bblid == ids[i]),] <- temp}}

int_sq <- scale(x[,23])^2
int_cub <- scale(x[,23])^3

#regressing out age and practice
ACC_apr <- scale(winsor(lm(pvrt_pc~Age+Age_Squared+Age_Cubed+TP+TP_Squared+TP_Cubed+int+int_sq+int_cub,data=x)$residuals,trim=0.005))
RT_apr <- scale(winsor(lm(pvrt_rtcr~Age+Age_Squared+Age_Cubed+TP+TP_Squared+TP_Cubed+int+int_sq+int_cub,data=x)$residuals,trim=0.005))

x <- data.frame(x,ACC_apr,RT_apr)

x1 <- x[which(x$timepoint == 1),]
x2 <- x[which(x$timepoint == 2),]
x3 <- x[which(x$timepoint == 3),]
x4 <- x[which(x$timepoint == 4),]
x5 <- x[which(x$timepoint == 5),]

x <- merge(x1,x2,by="bblid",all=TRUE)
x <- merge(x,x3,by="bblid",all=TRUE)
x <- merge(x,x4,by="bblid",all=TRUE)
x <- merge(x,x5,by="bblid",all=TRUE)

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
"pvrt_genus_TIME_1",
"pvrt_valid_TIME_1",
"pvrt_cr_TIME_1",
"pvrt_pc_TIME_1",
"pvrt_rtcr_TIME_1",
"timepoint_TIME_1",
"ntimepoints_TIME_1",
"ACC_ar_TIME_1",
"RT_ar_TIME_1",
"int_TIME_1",
"ACC_apr_TIME_1",
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
"pvrt_genus_TIME_2",
"pvrt_valid_TIME_2",
"pvrt_cr_TIME_2",
"pvrt_pc_TIME_2",
"pvrt_rtcr_TIME_2",
"timepoint_TIME_2",
"ntimepoints_TIME_2",
"ACC_ar_TIME_2",
"RT_ar_TIME_2",
"int_TIME_2",
"ACC_apr_TIME_2",
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
"pvrt_genus_TIME_3",
"pvrt_valid_TIME_3",
"pvrt_cr_TIME_3",
"pvrt_pc_TIME_3",
"pvrt_rtcr_TIME_3",
"timepoint_TIME_3",
"ntimepoints_TIME_3",
"ACC_ar_TIME_3",
"RT_ar_TIME_3",
"int_TIME_3",
"ACC_apr_TIME_3",
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
"pvrt_genus_TIME_4",
"pvrt_valid_TIME_4",
"pvrt_cr_TIME_4",
"pvrt_pc_TIME_4",
"pvrt_rtcr_TIME_4",
"timepoint_TIME_4",
"ntimepoints_TIME_4",
"ACC_ar_TIME_4",
"RT_ar_TIME_4",
"int_TIME_4",
"ACC_apr_TIME_4",
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
"pvrt_genus_TIME_5",
"pvrt_valid_TIME_5",
"pvrt_cr_TIME_5",
"pvrt_pc_TIME_5",
"pvrt_rtcr_TIME_5",
"timepoint_TIME_5",
"ntimepoints_TIME_5",
"ACC_ar_TIME_5",
"RT_ar_TIME_5",
"int_TIME_5",
"ACC_apr_TIME_5",
"RT_apr_TIME_5")

ACC <- data.frame(
x$pvrt_pc_TIME_1,
x$pvrt_pc_TIME_2,
x$pvrt_pc_TIME_3,
x$pvrt_pc_TIME_4,
x$pvrt_pc_TIME_5)

RT <- data.frame(
x$pvrt_rtcr_TIME_1,
x$pvrt_rtcr_TIME_2,
x$pvrt_rtcr_TIME_3,
x$pvrt_rtcr_TIME_4,
x$pvrt_rtcr_TIME_5)

ACC_ar <- data.frame(
x$ACC_ar_TIME_1,
x$ACC_ar_TIME_2,
x$ACC_ar_TIME_3,
x$ACC_ar_TIME_4,
x$ACC_ar_TIME_5)

RT_ar <- data.frame(
x$RT_ar_TIME_1,
x$RT_ar_TIME_2,
x$RT_ar_TIME_3,
x$RT_ar_TIME_4,
x$RT_ar_TIME_5)

ACC_apr <- data.frame(
x$ACC_apr_TIME_1,
x$ACC_apr_TIME_2,
x$ACC_apr_TIME_3,
x$ACC_apr_TIME_4,
x$ACC_apr_TIME_5)

RT_apr <- data.frame(
x$RT_apr_TIME_1,
x$RT_apr_TIME_2,
x$RT_apr_TIME_3,
x$RT_apr_TIME_4,
x$RT_apr_TIME_5)

res <- matrix(NA,4,4)

for (i in 2:5) {
res[1,(i-1)] <- icc(ACC_apr[,1:i],type="agreement",model="twoway")$value
res[2,(i-1)] <- icc(RT_apr[,1:i],type="agreement",model="twoway")$value
res[3,(i-1)] <- icc(ACC_ar[,1:i],type="agreement",model="twoway")$value
res[4,(i-1)] <- icc(RT_ar[,1:i],type="agreement",model="twoway")$value
} 

header <- c("2 Timepoints","3 Timepoints","4 Timepoints","5 Timepoints")
TestRetest_for_PVRT <- c("","Accuracy (age- & practice-regressed)","Response Time (age- & practice-regressed)","Accuracy (age-regressed)","Response Time (age-regressed)")

pdf("TestRetest_PVRT.pdf",width=15,height=15)
grid.table(cbind(TestRetest_for_PVRT,rbind(header,round(res,3))))
pairs.panels(ACC_apr,lm=TRUE)
pairs.panels(RT_apr,lm=TRUE)
dev.off()




