
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

##############################################################################
# CPF
##############################################################################

x <- read.csv("CNB Longitudinal/cpf_20191113.csv")

x[,4] <- as.numeric(ymd(x[,4])) # date of test
colnames(x)[4] <- "Time"
colnames(x)[3] <- "bblid"
colnames(x)[11] <- "Sex"
x <- x[order(x$bblid, x$timepoint),]

ids <- unique(x$bblid)
for (i in 1:length(ids)) {
temp <- x[which(x$bblid == ids[i]),]
if (is.element("N",temp[,6]) | is.element("N",temp[,15])) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])}
#try(if (min(temp$cpf_cr,na.rm=TRUE) < 10) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])})
try(if (max(temp$cpf_w_rtcr,na.rm=TRUE) > 10000) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])})
x[which(x$bblid == ids[i]),] <- temp}  

x <- x[which(x$ntimepoints > 1),]   # best is > 1
x <- x[which(x$timepoint < 5),]     # best is < 6
x <- x[which(is.na(x$cnbAgemonths) == FALSE),]

ids <- unique(x$bblid)
for (i in 1:length(ids)) {
temp <- x[which(x$bblid == ids[i]),]
temp[,18] <- seq(1:nrow(temp))
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

raw_raw <- x[,16:17]
x[,16:17] <- scale(x[,16:17])

ACC_r <- scale(winsor(lm(cpf_cr~Age+Age_Squared+Age_Cubed,data=x)$residuals,trim=0.005))
RT_r <- scale(winsor(lm(cpf_w_rtcr~Age+Age_Squared+Age_Cubed,data=x)$residuals,trim=0.005))

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
if (temp[j,18] == 4) {try(temp[j,22] <- temp[j,4] - temp[which(temp$timepoint == 3),4])}
if (temp[j,18] == 3) {try(temp[j,22] <- temp[j,4] - temp[which(temp$timepoint == 2),4])}
if (temp[j,18] == 2) {try(temp[j,22] <- temp[j,4])}
if (temp[j,18] == 1) {try(temp[j,22] <- NA)}
x[which(x$bblid == ids[i]),] <- temp}}

max_int <- max(x$int,na.rm=TRUE)

int_w_t1 <- x$int
int_w_t1[is.na(int_w_t1) == TRUE] <- max_int
x <- data.frame(x,int_w_t1)

int_sq <- scale(x[,23])^2
int_cub <- scale(x[,23])^3

#regressing out age and practice
ACC_apr <- scale(winsor(lm(cpf_cr~Age+Age_Squared+Age_Cubed+TP+TP_Squared+TP_Cubed+int_w_t1+int_sq+int_cub,data=x)$residuals,trim=0.005))
RT_apr <- scale(winsor(lm(cpf_w_rtcr~Age+Age_Squared+Age_Cubed+TP+TP_Squared+TP_Cubed+int_w_t1+int_sq+int_cub,data=x)$residuals,trim=0.005))

x$cnbAgemonths <- x$cnbAgemonths/12
x$Time <- x$Time/365.25
Time_Squared <- scale(x$Time)^2
Age <- x$cnbAgemonths
int_ord <- as.factor(quantcut(x$int,3))

x <- data.frame(x[,c(1,3:4,10:11,14:18,20:23)],ACC_apr,RT_apr,Age,"CPF",raw_raw)

colnames(x) <- c("	datasetid	","	bblid	","	Time	","	meducation	","	Sex	","	genus	",
"	valid_code	","	ACC_raw	","	RT_raw	","	timepoint	","	ACC_ar	","	RT_ar	","	int	",
"	int_w_t1	","	ACC_apr	","	RT_apr	","	Age	","Test","raw_original_acc","raw_original_RT")

cpf <- x


##############################################################################
# ER40
##############################################################################

x <- read.csv("CNB Longitudinal/er40_20191113.csv")

x[,4] <- as.numeric(ymd(x[,4])) # date of test
colnames(x)[4] <- "Time"
colnames(x)[3] <- "bblid"
colnames(x)[11] <- "Sex"
x <- x[order(x$bblid, x$timepoint),]

ids <- unique(x$bblid)
for (i in 1:length(ids)) {
temp <- x[which(x$bblid == ids[i]),]
if (is.element("N",temp[,6]) | is.element("N",temp[,15])) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])}
try(if (min(temp$er40_cr,na.rm=TRUE) < 10) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])})
try(if (max(temp$er40_rtcr,na.rm=TRUE) > 6000) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])})
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

#set.seed(2)
#temp <- amelia(x[,16:19], m=1)$imputations[[1]]
#x[,16:19] <- temp

raw_raw <- x[,16:17]
x[,16:17] <- scale(x[,16:17])

ACC_r <- scale(winsor(lm(er40_cr~Age+Age_Squared+Age_Cubed,data=x)$residuals,trim=0.005))
RT_r <- scale(winsor(lm(er40_rtcr~Age+Age_Squared+Age_Cubed,data=x)$residuals,trim=0.005))

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
if (temp[j,18] == 6) {try(temp[j,22] <- temp[j,4] - temp[which(temp$timepoint == 5),4])}
if (temp[j,18] == 5) {try(temp[j,22] <- temp[j,4] - temp[which(temp$timepoint == 4),4])}
if (temp[j,18] == 4) {try(temp[j,22] <- temp[j,4] - temp[which(temp$timepoint == 3),4])}
if (temp[j,18] == 3) {try(temp[j,22] <- temp[j,4] - temp[which(temp$timepoint == 2),4])}
if (temp[j,18] == 2) {try(temp[j,22] <- temp[j,4])}
if (temp[j,18] == 1) {try(temp[j,22] <- NA)}
x[which(x$bblid == ids[i]),] <- temp}}

max_int <- max(x$int,na.rm=TRUE)

int_w_t1 <- x$int
int_w_t1[is.na(int_w_t1) == TRUE] <- max_int
x <- data.frame(x,int_w_t1)

int_sq <- scale(x[,23])^2
int_cub <- scale(x[,23])^3

#regressing out age and practice
ACC_apr <- scale(winsor(lm(er40_cr~Age+Age_Squared+Age_Cubed+TP+TP_Squared+TP_Cubed+int_w_t1+int_sq+int_cub,data=x)$residuals,trim=0.005))
RT_apr <- scale(winsor(lm(er40_rtcr~Age+Age_Squared+Age_Cubed+TP+TP_Squared+TP_Cubed+int_w_t1+int_sq+int_cub,data=x)$residuals,trim=0.005))

x$cnbAgemonths <- x$cnbAgemonths/12
x$Time <- x$Time/365.25
Time_Squared <- scale(x$Time)^2
Age <- x$cnbAgemonths
int_ord <- as.factor(quantcut(x$int,3))

x <- data.frame(x[,c(1,3:4,10:11,14:18,20:23)],ACC_apr,RT_apr,Age,"ER40",raw_raw)

colnames(x) <- c("	datasetid	","	bblid	","	Time	","	meducation	","	Sex	","	genus	",
"	valid_code	","	ACC_raw	","	RT_raw	","	timepoint	","	ACC_ar	","	RT_ar	","	int	",
"	int_w_t1	","	ACC_apr	","	RT_apr	","	Age	","Test","raw_original_acc","raw_original_RT")

er40 <- x


##############################################################################
# PVRT
##############################################################################

x <- read.csv("CNB Longitudinal/pvrt_20191113.csv")

x[,4] <- as.numeric(ymd(x[,4])) # date of test
colnames(x)[4] <- "Time"
colnames(x)[3] <- "bblid"
colnames(x)[11] <- "Sex"
x <- x[order(x$bblid, x$timepoint),]

x <- x[grepl("KSPVRT",x$pvrt_genus),]

ids <- unique(x$bblid)
for (i in 1:length(ids)) {
temp <- x[which(x$bblid == ids[i]),]
if (is.element("N",temp[,6]) | is.element("N",temp[,15])) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])}
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

raw_raw <- x[,17:18]
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

int_w_t1 <- x$int
int_w_t1[is.na(int_w_t1) == TRUE] <- max_int
x <- data.frame(x,int_w_t1)

int_sq <- scale(x[,24])^2
int_cub <- scale(x[,24])^3

#regressing out age and practice
ACC_apr <- scale(winsor(lm(pvrt_pc~Age+Age_Squared+Age_Cubed+TP+TP_Squared+TP_Cubed+int_w_t1+int_sq+int_cub,data=x)$residuals,trim=0.005))
RT_apr <- scale(winsor(lm(pvrt_rtcr~Age+Age_Squared+Age_Cubed+TP+TP_Squared+TP_Cubed+int_w_t1+int_sq+int_cub,data=x)$residuals,trim=0.005))

x$cnbAgemonths <- x$cnbAgemonths/12
x$Time <- x$Time/365.25
Time_Squared <- scale(x$Time)^2
Age <- x$cnbAgemonths
int_ord <- as.factor(quantcut(x$int,3))

x <- data.frame(x[,c(1,3:4,10:11,14:15,17:19,21:24)],ACC_apr,RT_apr,Age,"PVRT",raw_raw)

colnames(x) <- c("	datasetid	","	bblid	","	Time	","	meducation	","	Sex	","	genus	",
"	valid_code	","	ACC_raw	","	RT_raw	","	timepoint	","	ACC_ar	","	RT_ar	","	int	",
"	int_w_t1	","	ACC_apr	","	RT_apr	","	Age	","Test","raw_original_acc","raw_original_RT")

pvrt <- x


##############################################################################
# CPT
##############################################################################

x <- read.csv("CNB Longitudinal/cpt_20191113.csv")

x[,4] <- as.numeric(ymd(x[,4])) # date of test
colnames(x)[4] <- "Time"
colnames(x)[3] <- "bblid"
colnames(x)[11] <- "Sex"
x <- x[order(x$bblid, x$timepoint),]

ids <- unique(x$bblid)
for (i in 1:length(ids)) {
temp <- x[which(x$bblid == ids[i]),]
if (is.element("N",temp[,6]) | is.element("N",temp[,15])) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])}
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

raw_raw <- data.frame(qnorm(x$cpt_ptp-0.001) - qnorm(x$cpt_pfp+0.001),x[,17])
dprime <- scale(qnorm(x$cpt_ptp-0.001) - qnorm(x$cpt_pfp+0.001))

x[,16:19] <- scale(x[,16:19])

PTP_r <- scale(winsor(lm(cpt_ptp~Age+Age_Squared+Age_Cubed,data=x)$residuals,trim=0.005))
dmod <- lm(dprime~Age+Age_Squared+Age_Cubed,data=x,na.action=na.exclude)
ACC_r <- scale(winsor(as.matrix(residuals(dmod,na.action=na.exclude)),trim=0.005))
TPRT_r <- scale(winsor(lm(cpt_tprt~Age+Age_Squared+Age_Cubed,data=x)$residuals,trim=0.005))
FPRT_r <- scale(winsor(lm(cpt_fprt~Age+Age_Squared+Age_Cubed,data=x)$residuals,trim=0.005))
RT_r <- scale(winsor(lm(cpt_tprt~Age+Age_Squared+Age_Cubed,data=x)$residuals,trim=0.005))

x <- data.frame(x,PTP_r,ACC_r,TPRT_r,FPRT_r,RT_r)

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
if (temp[j,20] == 6) {try(temp[j,27] <- temp[j,4] - temp[which(temp$timepoint == 5),4])}
if (temp[j,20] == 5) {try(temp[j,27] <- temp[j,4] - temp[which(temp$timepoint == 4),4])}
if (temp[j,20] == 4) {try(temp[j,27] <- temp[j,4] - temp[which(temp$timepoint == 3),4])}
if (temp[j,20] == 3) {try(temp[j,27] <- temp[j,4] - temp[which(temp$timepoint == 2),4])}
if (temp[j,20] == 2) {try(temp[j,27] <- temp[j,4])}
if (temp[j,20] == 1) {try(temp[j,27] <- NA)}
x[which(x$bblid == ids[i]),] <- temp}}

max_int <- max(x$int,na.rm=TRUE)

int_w_t1 <- x$int
int_w_t1[is.na(int_w_t1) == TRUE] <- max_int
x <- data.frame(x,int_w_t1)

int_sq <- scale(x[,28])^2
int_cub <- scale(x[,28])^3

#regressing out age and practice
ACC_apr <- scale(winsor(as.matrix(lm(dprime~Age+Age_Squared+Age_Cubed+TP+TP_Squared+TP_Cubed+int_w_t1+int_sq+int_cub,data=x)$residuals),trim=0.005))
RT_apr <- scale(winsor(lm(cpt_tprt~Age+Age_Squared+Age_Cubed+TP+TP_Squared+TP_Cubed+int_w_t1+int_sq+int_cub,data=x)$residuals,trim=0.005))

x$cnbAgemonths <- x$cnbAgemonths/12
x$Time <- x$Time/365.25
Time_Squared <- scale(x$Time)^2
Age <- x$cnbAgemonths
int_ord <- as.factor(quantcut(x$int,3))

x$cpt_ptp <- dprime

x <- data.frame(x[,c(1,3:4,10:11,14:17,20,23,26:28)],ACC_apr,RT_apr,Age,"CPT",raw_raw)

colnames(x) <- c("	datasetid	","	bblid	","	Time	","	meducation	","	Sex	","	genus	",
"	valid_code	","	ACC_raw	","	RT_raw	","	timepoint	","	ACC_ar	","	RT_ar	","	int	",
"	int_w_t1	","	ACC_apr	","	RT_apr	","	Age	","Test","raw_original_acc","raw_original_RT")

cpt <- x


##############################################################################
# PCET
##############################################################################

x <- read.csv("CNB Longitudinal/pcet_20191113.csv")

dates <- as.numeric(ymd(x[,4])) # date of test
if (sum(!is.na(dates)) < 10) {dates <- as.numeric(mdy(x[,4]))}
x[,4] <- dates
colnames(x)[4] <- "Time"
colnames(x)[3] <- "bblid"
colnames(x)[11] <- "Sex"
x <- x[order(x$bblid, x$timepoint),]

ids <- unique(x$bblid)
for (i in 1:length(ids)) {
temp <- x[which(x$bblid == ids[i]),]
if (is.element("N",temp[,6]) | is.element("N",temp[,15])) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])}
#try(if (min(temp$pcet_acc2,na.rm=TRUE) < 10) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])})
try(if (max(temp$pcet_rtcr,na.rm=TRUE) > 15000) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])})
x[which(x$bblid == ids[i]),] <- temp}  

x <- x[which(x$ntimepoints > 1),]   # best is > 1
x <- x[which(x$timepoint < 5),]     # best is < 6
x <- x[which(is.na(x$cnbAgemonths) == FALSE),]

ids <- unique(x$bblid)
for (i in 1:length(ids)) {
temp <- x[which(x$bblid == ids[i]),]
temp[,18] <- seq(1:nrow(temp))
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

raw_raw <- x[,16:17]
x[,16:17] <- scale(x[,16:17])

ACC_r <- scale(winsor(lm(pcet_acc2~Age+Age_Squared+Age_Cubed,data=x)$residuals,trim=0.005))
RT_r <- scale(winsor(lm(pcet_rtcr~Age+Age_Squared+Age_Cubed,data=x)$residuals,trim=0.005))

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
if (temp[j,18] == 4) {try(temp[j,22] <- temp[j,4] - temp[which(temp$timepoint == 3),4])}
if (temp[j,18] == 3) {try(temp[j,22] <- temp[j,4] - temp[which(temp$timepoint == 2),4])}
if (temp[j,18] == 2) {try(temp[j,22] <- temp[j,4])}
if (temp[j,18] == 1) {try(temp[j,22] <- NA)}
x[which(x$bblid == ids[i]),] <- temp}}

max_int <- max(x$int,na.rm=TRUE)

int_w_t1 <- x$int
int_w_t1[is.na(int_w_t1) == TRUE] <- max_int
x <- data.frame(x,int_w_t1)

int_sq <- scale(x[,23])^2
int_cub <- scale(x[,23])^3

#regressing out age and practice
ACC_apr <- scale(winsor(lm(pcet_acc2~Age+Age_Squared+Age_Cubed+TP+TP_Squared+TP_Cubed+int_w_t1+int_sq+int_cub,data=x)$residuals,trim=0.005))
RT_apr <- scale(winsor(lm(pcet_rtcr~Age+Age_Squared+Age_Cubed+TP+TP_Squared+TP_Cubed+int_w_t1+int_sq+int_cub,data=x)$residuals,trim=0.005))

x$cnbAgemonths <- x$cnbAgemonths/12
x$Time <- x$Time/365.25
Time_Squared <- scale(x$Time)^2
Age <- x$cnbAgemonths
int_ord <- as.factor(quantcut(x$int,3))

x <- data.frame(x[,c(1,3:4,10:11,14:18,20:23)],ACC_apr,RT_apr,Age,"PCET",raw_raw)

colnames(x) <- c("	datasetid	","	bblid	","	Time	","	meducation	","	Sex	","	genus	",
"	valid_code	","	ACC_raw	","	RT_raw	","	timepoint	","	ACC_ar	","	RT_ar	","	int	",
"	int_w_t1	","	ACC_apr	","	RT_apr	","	Age	","Test","raw_original_acc","raw_original_RT")

pcet <- x


##############################################################################
# PLOT
##############################################################################

x <- read.csv("CNB Longitudinal/plot_20191113.csv")

dates <- as.numeric(ymd(x[,4])) # date of test
if (sum(!is.na(dates)) < 10) {dates <- as.numeric(mdy(x[,4]))}
x[,4] <- dates
colnames(x)[4] <- "Time"
colnames(x)[3] <- "bblid"
colnames(x)[11] <- "Sex"
x <- x[order(x$bblid, x$timepoint),]

ids <- unique(x$bblid)
for (i in 1:length(ids)) {
temp <- x[which(x$bblid == ids[i]),]
if (is.element("N",temp[,6]) | is.element("N",temp[,15])) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])}
#try(if (min(temp$plot_pc,na.rm=TRUE) < 10) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])})
try(if (max(temp$plot_rtcr,na.rm=TRUE) > 30000) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])})
x[which(x$bblid == ids[i]),] <- temp}  

x <- x[which(x$ntimepoints > 1),]   # best is > 1
x <- x[which(x$timepoint < 5),]     # best is < 6
x <- x[which(is.na(x$cnbAgemonths) == FALSE),]
x <- x[which(is.na(x$plot_pc) == FALSE),]

ids <- unique(x$bblid)
for (i in 1:length(ids)) {
temp <- x[which(x$bblid == ids[i]),]
temp[,18] <- seq(1:nrow(temp))
x[which(x$bblid == ids[i]),] <- temp}

Age <- scale(x$cnbAgemonths)
Age_Squared <- Age^2
Age_Cubed <- Age^3
TP <- scale(x$timepoint)
TP_Squared <- scale(x$timepoint)^2
TP_Cubed <- scale(x$timepoint)^3

set.seed(2)
temp <- amelia(x[,16:17], m=1)$imputations[[1]]
x[,16:17] <- temp

raw_raw <- x[,16:17]
x[,16:17] <- scale(x[,16:17])

ACC_r <- scale(winsor(lm(plot_pc~Age+Age_Squared+Age_Cubed,data=x)$residuals,trim=0.005))
RT_r <- scale(winsor(lm(plot_rtcr~Age+Age_Squared+Age_Cubed,data=x)$residuals,trim=0.005))

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
if (temp[j,18] == 4) {try(temp[j,22] <- temp[j,4] - temp[which(temp$timepoint == 3),4])}
if (temp[j,18] == 3) {try(temp[j,22] <- temp[j,4] - temp[which(temp$timepoint == 2),4])}
if (temp[j,18] == 2) {try(temp[j,22] <- temp[j,4])}
if (temp[j,18] == 1) {try(temp[j,22] <- NA)}
x[which(x$bblid == ids[i]),] <- temp}}

max_int <- max(x$int,na.rm=TRUE)

int_w_t1 <- x$int
int_w_t1[is.na(int_w_t1) == TRUE] <- max_int
x <- data.frame(x,int_w_t1)

int_sq <- scale(x[,23])^2
int_cub <- scale(x[,23])^3

#regressing out age and practice
ACC_apr <- scale(winsor(lm(plot_pc~Age+Age_Squared+Age_Cubed+TP+TP_Squared+TP_Cubed+int_w_t1+int_sq+int_cub,data=x)$residuals,trim=0.005))
RT_apr <- scale(winsor(lm(plot_rtcr~Age+Age_Squared+Age_Cubed+TP+TP_Squared+TP_Cubed+int_w_t1+int_sq+int_cub,data=x)$residuals,trim=0.005))

x$cnbAgemonths <- x$cnbAgemonths/12
x$Time <- x$Time/365.25
Time_Squared <- scale(x$Time)^2
Age <- x$cnbAgemonths
int_ord <- as.factor(quantcut(x$int,3))

x <- data.frame(x[,c(1,3:4,10:11,14:18,20:23)],ACC_apr,RT_apr,Age,"PLOT",raw_raw)

colnames(x) <- c("	datasetid	","	bblid	","	Time	","	meducation	","	Sex	","	genus	",
"	valid_code	","	ACC_raw	","	RT_raw	","	timepoint	","	ACC_ar	","	RT_ar	","	int	",
"	int_w_t1	","	ACC_apr	","	RT_apr	","	Age	","Test","raw_original_acc","raw_original_RT")

plOt <- x


##############################################################################
# PMAT
##############################################################################

x <- read.csv("CNB Longitudinal/pmat_20191113.csv")

x[,4] <- as.numeric(ymd(x[,4])) # date of test
colnames(x)[4] <- "Time"
colnames(x)[3] <- "bblid"
colnames(x)[11] <- "Sex"
x <- x[order(x$bblid, x$timepoint),]

ids <- unique(x$bblid)
for (i in 1:length(ids)) {
temp <- x[which(x$bblid == ids[i]),]
if (is.element("N",temp[,6]) | is.element("N",temp[,15])) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])}
#try(if (min(temp$pmat_pc,na.rm=TRUE) < 10) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])})
try(if (max(temp$pmat_rtcr,na.rm=TRUE) > 35000) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])})
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

set.seed(2)
temp <- amelia(x[,16:18], m=1)$imputations[[1]]
x[,16:18] <- temp

raw_raw <- x[,c(16,18)]
x[,16:18] <- scale(x[,16:18])

ACC_r <- scale(winsor(lm(pmat_pc~Age+Age_Squared+Age_Cubed,data=x)$residuals,trim=0.005))
RT_r <- scale(winsor(lm(pmat_rtcr~Age+Age_Squared+Age_Cubed,data=x)$residuals,trim=0.005))

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

int_w_t1 <- x$int
int_w_t1[is.na(int_w_t1) == TRUE] <- max_int
x <- data.frame(x,int_w_t1)

int_sq <- scale(x[,24])^2
int_cub <- scale(x[,24])^3

#regressing out age and practice
ACC_apr <- scale(winsor(lm(pmat_pc~Age+Age_Squared+Age_Cubed+TP+TP_Squared+TP_Cubed+int_w_t1+int_sq+int_cub,data=x)$residuals,trim=0.005))
RT_apr <- scale(winsor(lm(pmat_rtcr~Age+Age_Squared+Age_Cubed+TP+TP_Squared+TP_Cubed+int_w_t1+int_sq+int_cub,data=x)$residuals,trim=0.005))

x$cnbAgemonths <- x$cnbAgemonths/12
x$Time <- x$Time/365.25
Time_Squared <- scale(x$Time)^2
Age <- x$cnbAgemonths
int_ord <- as.factor(quantcut(x$int,3))

x <- data.frame(x[,c(1,3:4,10:11,14:16,18:19,21:24)],ACC_apr,RT_apr,Age,"PMAT",raw_raw)

colnames(x) <- c("	datasetid	","	bblid	","	Time	","	meducation	","	Sex	","	genus	",
"	valid_code	","	ACC_raw	","	RT_raw	","	timepoint	","	ACC_ar	","	RT_ar	","	int	",
"	int_w_t1	","	ACC_apr	","	RT_apr	","	Age	","Test","raw_original_acc","raw_original_RT")

pmat <- x


##############################################################################
# VOLT
##############################################################################

x <- read.csv("CNB Longitudinal/volt_20191113.csv")

x[,4] <- as.numeric(ymd(x[,4])) # date of test
colnames(x)[4] <- "Time"
colnames(x)[3] <- "bblid"
colnames(x)[11] <- "Sex"
x <- x[order(x$bblid, x$timepoint),]

ids <- unique(x$bblid)
for (i in 1:length(ids)) {
temp <- x[which(x$bblid == ids[i]),]
if (is.element("N",temp[,6]) | is.element("N",temp[,15])) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])}
#try(if (min(temp$volt_cr,na.rm=TRUE) < 10) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])})
try(if (max(temp$volt_w_rtcr,na.rm=TRUE) > 6000) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])})
x[which(x$bblid == ids[i]),] <- temp}  

x <- x[which(x$ntimepoints > 1),]   # best is > 1
x <- x[which(x$timepoint < 5),]     # best is < 6
x <- x[which(is.na(x$cnbAgemonths) == FALSE),]

ids <- unique(x$bblid)
for (i in 1:length(ids)) {
temp <- x[which(x$bblid == ids[i]),]
temp[,18] <- seq(1:nrow(temp))
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

raw_raw <- x[,16:17]
x[,16:17] <- scale(x[,16:17])

ACC_r <- scale(winsor(lm(volt_cr~Age+Age_Squared+Age_Cubed,data=x)$residuals,trim=0.005))
RT_r <- scale(winsor(lm(volt_w_rtcr~Age+Age_Squared+Age_Cubed,data=x)$residuals,trim=0.005))

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
if (temp[j,18] == 4) {try(temp[j,22] <- temp[j,4] - temp[which(temp$timepoint == 3),4])}
if (temp[j,18] == 3) {try(temp[j,22] <- temp[j,4] - temp[which(temp$timepoint == 2),4])}
if (temp[j,18] == 2) {try(temp[j,22] <- temp[j,4])}
if (temp[j,18] == 1) {try(temp[j,22] <- NA)}
x[which(x$bblid == ids[i]),] <- temp}}

max_int <- max(x$int,na.rm=TRUE)

int_w_t1 <- x$int
int_w_t1[is.na(int_w_t1) == TRUE] <- max_int
x <- data.frame(x,int_w_t1)

int_sq <- scale(x[,23])^2
int_cub <- scale(x[,23])^3

#regressing out age and practice
ACC_apr <- scale(winsor(lm(volt_cr~Age+Age_Squared+Age_Cubed+TP+TP_Squared+TP_Cubed+int_w_t1+int_sq+int_cub,data=x)$residuals,trim=0.005))
RT_apr <- scale(winsor(lm(volt_w_rtcr~Age+Age_Squared+Age_Cubed+TP+TP_Squared+TP_Cubed+int_w_t1+int_sq+int_cub,data=x)$residuals,trim=0.005))

x$cnbAgemonths <- x$cnbAgemonths/12
x$Time <- x$Time/365.25
Time_Squared <- scale(x$Time)^2
Age <- x$cnbAgemonths
int_ord <- as.factor(quantcut(x$int,3))

x <- data.frame(x[,c(1,3:4,10:11,14:18,20:23)],ACC_apr,RT_apr,Age,"VOLT",raw_raw)

colnames(x) <- c("	datasetid	","	bblid	","	Time	","	meducation	","	Sex	","	genus	",
"	valid_code	","	ACC_raw	","	RT_raw	","	timepoint	","	ACC_ar	","	RT_ar	","	int	",
"	int_w_t1	","	ACC_apr	","	RT_apr	","	Age	","Test","raw_original_acc","raw_original_RT")

volt <- x


##############################################################################
# NBACK
##############################################################################

x <- read.csv("CNB Longitudinal/lnb_20191113.csv")

x[,4] <- as.numeric(ymd(x[,4])) # date of test
colnames(x)[4] <- "Time"
colnames(x)[3] <- "bblid"
colnames(x)[11] <- "Sex"
x <- x[order(x$bblid, x$timepoint),]

ids <- unique(x$bblid)
for (i in 1:length(ids)) {
temp <- x[which(x$bblid == ids[i]),]
if (is.element("N",temp[,6]) | is.element("N",temp[,15]) | is.element("lnb2_135",temp[,14]) | is.element("lnb2_141",temp[,14])) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])}
try(if (min(temp$lnb_mpc,na.rm=TRUE) < 20) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])})
try(if (max(temp$lnb_mrtc,na.rm=TRUE) > 1500) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])})
x[which(x$bblid == ids[i]),] <- temp}  

x <- x[which(x$ntimepoints > 1),]   # best is > 1
x <- x[which(x$timepoint < 5),]     # best is < 6
x <- x[which(is.na(x$cnbAgemonths) == FALSE),]

ids <- unique(x$bblid)
for (i in 1:length(ids)) {
temp <- x[which(x$bblid == ids[i]),]
temp[,21] <- seq(1:nrow(temp))
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

raw_raw <- x[,17:18]
x[,16:20] <- scale(x[,16:20])

ACC_r <- scale(winsor(lm(lnb_mpc~Age+Age_Squared+Age_Cubed,data=x)$residuals,trim=0.005))
RT_r <- scale(winsor(lm(lnb_mrtc~Age+Age_Squared+Age_Cubed,data=x)$residuals,trim=0.005))

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
if (temp[j,21] == 4) {try(temp[j,25] <- temp[j,4] - temp[which(temp$timepoint == 3),4])}
if (temp[j,21] == 3) {try(temp[j,25] <- temp[j,4] - temp[which(temp$timepoint == 2),4])}
if (temp[j,21] == 2) {try(temp[j,25] <- temp[j,4])}
if (temp[j,21] == 1) {try(temp[j,25] <- NA)}
x[which(x$bblid == ids[i]),] <- temp}}

max_int <- max(x$int,na.rm=TRUE)

int_w_t1 <- x$int
int_w_t1[is.na(int_w_t1) == TRUE] <- max_int
x <- data.frame(x,int_w_t1)

int_sq <- scale(x[,26])^2
int_cub <- scale(x[,26])^3

#regressing out age and practice
ACC_apr <- scale(winsor(lm(lnb_mpc~Age+Age_Squared+Age_Cubed+TP+TP_Squared+TP_Cubed+int_w_t1+int_sq+int_cub,data=x)$residuals,trim=0.005))
RT_apr <- scale(winsor(lm(lnb_mrtc~Age+Age_Squared+Age_Cubed+TP+TP_Squared+TP_Cubed+int_w_t1+int_sq+int_cub,data=x)$residuals,trim=0.005))

x$cnbAgemonths <- x$cnbAgemonths/12
x$Time <- x$Time/365.25
Time_Squared <- scale(x$Time)^2
Age <- x$cnbAgemonths
int_ord <- as.factor(quantcut(x$int,3))

x <- data.frame(x[,c(1,3:4,10:11,14:15,17:18,21,23:26)],ACC_apr,RT_apr,Age,"NBACK",raw_raw)

colnames(x) <- c("	datasetid	","	bblid	","	Time	","	meducation	","	Sex	","	genus	",
"	valid_code	","	ACC_raw	","	RT_raw	","	timepoint	","	ACC_ar	","	RT_ar	","	int	",
"	int_w_t1	","	ACC_apr	","	RT_apr	","	Age	","Test","raw_original_acc","raw_original_RT")

nback <- x


##############################################################################
# MEDF
##############################################################################

x <- read.csv("CNB Longitudinal/medf_20191113.csv")

x[,4] <- as.numeric(ymd(x[,4])) # date of test
colnames(x)[4] <- "Time"
colnames(x)[3] <- "bblid"
colnames(x)[11] <- "Sex"
x <- x[order(x$bblid, x$timepoint),]

ids <- unique(x$bblid)
for (i in 1:length(ids)) {
temp <- x[which(x$bblid == ids[i]),]
if (is.element("N",temp[,6]) | is.element("N",temp[,15])) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])}
try(if (min(temp$medf_pc,na.rm=TRUE) < 20) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])})
try(if (max(temp$medf_rtcr,na.rm=TRUE) > 9000) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])})
x[which(x$bblid == ids[i]),] <- temp}  

x <- x[which(x$ntimepoints > 1),]   # best is > 1
x <- x[which(x$timepoint < 5),]     # best is < 6
x <- x[which(is.na(x$cnbAgemonths) == FALSE),]

ids <- unique(x$bblid)
for (i in 1:length(ids)) {
temp <- x[which(x$bblid == ids[i]),]
temp[,18] <- seq(1:nrow(temp))
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

raw_raw <- x[,16:17]
x[,16:17] <- scale(x[,16:17])

ACC_r <- scale(winsor(lm(medf_pc~Age+Age_Squared+Age_Cubed,data=x)$residuals,trim=0.005))
RT_r <- scale(winsor(lm(medf_rtcr~Age+Age_Squared+Age_Cubed,data=x)$residuals,trim=0.005))

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
if (temp[j,18] == 4) {try(temp[j,22] <- temp[j,4] - temp[which(temp$timepoint == 3),4])}
if (temp[j,18] == 3) {try(temp[j,22] <- temp[j,4] - temp[which(temp$timepoint == 2),4])}
if (temp[j,18] == 2) {try(temp[j,22] <- temp[j,4])}
if (temp[j,18] == 1) {try(temp[j,22] <- NA)}
x[which(x$bblid == ids[i]),] <- temp}}

max_int <- max(x$int,na.rm=TRUE)

int_w_t1 <- x$int
int_w_t1[is.na(int_w_t1) == TRUE] <- max_int
x <- data.frame(x,int_w_t1)

int_sq <- scale(x[,23])^2
int_cub <- scale(x[,23])^3

#regressing out age and practice
ACC_apr <- scale(winsor(lm(medf_pc~Age+Age_Squared+Age_Cubed+TP+TP_Squared+TP_Cubed+int_w_t1+int_sq+int_cub,data=x)$residuals,trim=0.005))
RT_apr <- scale(winsor(lm(medf_rtcr~Age+Age_Squared+Age_Cubed+TP+TP_Squared+TP_Cubed+int_w_t1+int_sq+int_cub,data=x)$residuals,trim=0.005))

x$cnbAgemonths <- x$cnbAgemonths/12
x$Time <- x$Time/365.25
Time_Squared <- scale(x$Time)^2
Age <- x$cnbAgemonths
int_ord <- as.factor(quantcut(x$int,3))

x <- data.frame(x[,c(1,3:4,10:11,14:18,20:23)],ACC_apr,RT_apr,Age,"MEDF",raw_raw)

colnames(x) <- c("	datasetid	","	bblid	","	Time	","	meducation	","	Sex	","	genus	",
"	valid_code	","	ACC_raw	","	RT_raw	","	timepoint	","	ACC_ar	","	RT_ar	","	int	",
"	int_w_t1	","	ACC_apr	","	RT_apr	","	Age	","Test","raw_original_acc","raw_original_RT")

medf <- x


##############################################################################
# ADT
##############################################################################

x <- read.csv("CNB Longitudinal/adt_20191113.csv")

x[,4] <- as.numeric(ymd(x[,4])) # date of test
colnames(x)[4] <- "Time"
colnames(x)[3] <- "bblid"
colnames(x)[11] <- "Sex"
x <- x[order(x$bblid, x$timepoint),]

ids <- unique(x$bblid)
for (i in 1:length(ids)) {
temp <- x[which(x$bblid == ids[i]),]
if (is.element("N",temp[,6]) | is.element("N",temp[,15])) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])}
try(if (min(temp$adt_pc,na.rm=TRUE) < 20) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])})
try(if (max(temp$adt_rtcr,na.rm=TRUE) > 7000) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])})
x[which(x$bblid == ids[i]),] <- temp}  

x <- x[which(x$ntimepoints > 1),]   # best is > 1
x <- x[which(x$timepoint < 5),]     # best is < 6
x <- x[which(is.na(x$cnbAgemonths) == FALSE),]

ids <- unique(x$bblid)
for (i in 1:length(ids)) {
temp <- x[which(x$bblid == ids[i]),]
temp[,18] <- seq(1:nrow(temp))
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

raw_raw <- x[,16:17]
x[,16:17] <- scale(x[,16:17])

ACC_r <- scale(winsor(lm(adt_pc~Age+Age_Squared+Age_Cubed,data=x)$residuals,trim=0.005))
RT_r <- scale(winsor(lm(adt_rtcr~Age+Age_Squared+Age_Cubed,data=x)$residuals,trim=0.005))

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
if (temp[j,18] == 4) {try(temp[j,22] <- temp[j,4] - temp[which(temp$timepoint == 3),4])}
if (temp[j,18] == 3) {try(temp[j,22] <- temp[j,4] - temp[which(temp$timepoint == 2),4])}
if (temp[j,18] == 2) {try(temp[j,22] <- temp[j,4])}
if (temp[j,18] == 1) {try(temp[j,22] <- NA)}
x[which(x$bblid == ids[i]),] <- temp}}

max_int <- max(x$int,na.rm=TRUE)

int_w_t1 <- x$int
int_w_t1[is.na(int_w_t1) == TRUE] <- max_int
x <- data.frame(x,int_w_t1)

int_sq <- scale(x[,23])^2
int_cub <- scale(x[,23])^3

#regressing out age and practice
ACC_apr <- scale(winsor(lm(adt_pc~Age+Age_Squared+Age_Cubed+TP+TP_Squared+TP_Cubed+int_w_t1+int_sq+int_cub,data=x)$residuals,trim=0.005))
RT_apr <- scale(winsor(lm(adt_rtcr~Age+Age_Squared+Age_Cubed+TP+TP_Squared+TP_Cubed+int_w_t1+int_sq+int_cub,data=x)$residuals,trim=0.005))

x$cnbAgemonths <- x$cnbAgemonths/12
x$Time <- x$Time/365.25
Time_Squared <- scale(x$Time)^2
Age <- x$cnbAgemonths
int_ord <- as.factor(quantcut(x$int,3))

x <- data.frame(x[,c(1,3:4,10:11,14:18,20:23)],ACC_apr,RT_apr,Age,"ADT",raw_raw)

colnames(x) <- c("	datasetid	","	bblid	","	Time	","	meducation	","	Sex	","	genus	",
"	valid_code	","	ACC_raw	","	RT_raw	","	timepoint	","	ACC_ar	","	RT_ar	","	int	",
"	int_w_t1	","	ACC_apr	","	RT_apr	","	Age	","Test","raw_original_acc","raw_original_RT")

adt <- x


##############################################################################
# CPW
##############################################################################

x <- read.csv("CNB Longitudinal/cpw_20191113.csv")

x[,4] <- as.numeric(ymd(x[,4])) # date of test
colnames(x)[4] <- "Time"
colnames(x)[3] <- "bblid"
colnames(x)[11] <- "Sex"
x <- x[order(x$bblid, x$timepoint),]

ids <- unique(x$bblid)
for (i in 1:length(ids)) {
temp <- x[which(x$bblid == ids[i]),]
if (is.element("N",temp[,6]) | is.element("N",temp[,15]) | is.element("cpw_a",temp[,14]) | is.element("cpw_b",temp[,14])) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])}
try(if (min(temp$cpw_cr,na.rm=TRUE) < 15) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])})
try(if (max(temp$cpw_w_rtcr,na.rm=TRUE) > 6000) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])})
x[which(x$bblid == ids[i]),] <- temp}  

x <- x[which(x$ntimepoints > 1),]   # best is > 1
x <- x[which(x$timepoint < 5),]     # best is < 6
x <- x[which(is.na(x$cnbAgemonths) == FALSE),]

ids <- unique(x$bblid)
for (i in 1:length(ids)) {
temp <- x[which(x$bblid == ids[i]),]
temp[,18] <- seq(1:nrow(temp))
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

raw_raw <- x[,16:17]
x[,16:17] <- scale(x[,16:17])

ACC_r <- scale(winsor(lm(cpw_cr~Age+Age_Squared+Age_Cubed,data=x)$residuals,trim=0.005))
RT_r <- scale(winsor(lm(cpw_w_rtcr~Age+Age_Squared+Age_Cubed,data=x)$residuals,trim=0.005))

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
if (temp[j,18] == 4) {try(temp[j,22] <- temp[j,4] - temp[which(temp$timepoint == 3),4])}
if (temp[j,18] == 3) {try(temp[j,22] <- temp[j,4] - temp[which(temp$timepoint == 2),4])}
if (temp[j,18] == 2) {try(temp[j,22] <- temp[j,4])}
if (temp[j,18] == 1) {try(temp[j,22] <- NA)}
x[which(x$bblid == ids[i]),] <- temp}}

max_int <- max(x$int,na.rm=TRUE)

int_w_t1 <- x$int
int_w_t1[is.na(int_w_t1) == TRUE] <- max_int
x <- data.frame(x,int_w_t1)

int_sq <- scale(x[,23])^2
int_cub <- scale(x[,23])^3

#regressing out age and practice
ACC_apr <- scale(winsor(lm(cpw_cr~Age+Age_Squared+Age_Cubed+TP+TP_Squared+TP_Cubed+int_w_t1+int_sq+int_cub,data=x)$residuals,trim=0.005))
RT_apr <- scale(winsor(lm(cpw_w_rtcr~Age+Age_Squared+Age_Cubed+TP+TP_Squared+TP_Cubed+int_w_t1+int_sq+int_cub,data=x)$residuals,trim=0.005))

x$cnbAgemonths <- x$cnbAgemonths/12
x$Time <- x$Time/365.25
Time_Squared <- scale(x$Time)^2
Age <- x$cnbAgemonths
int_ord <- as.factor(quantcut(x$int,3))

x <- data.frame(x[,c(1,3:4,10:11,14:18,20:23)],ACC_apr,RT_apr,Age,"CPW",raw_raw)

colnames(x) <- c("	datasetid	","	bblid	","	Time	","	meducation	","	Sex	","	genus	",
"	valid_code	","	ACC_raw	","	RT_raw	","	timepoint	","	ACC_ar	","	RT_ar	","	int	",
"	int_w_t1	","	ACC_apr	","	RT_apr	","	Age	","Test","raw_original_acc","raw_original_RT")

cpw <- x


##############################################################################
# MPRAXIS
##############################################################################

x <- read.csv("CNB Longitudinal/mpraxis_20191113.csv")

x[,4] <- as.numeric(ymd(x[,4])) # date of test
colnames(x)[4] <- "Time"
colnames(x)[3] <- "bblid"
colnames(x)[11] <- "Sex"
x <- x[order(x$bblid, x$timepoint),]

ids <- unique(x$bblid)
for (i in 1:length(ids)) {
temp <- x[which(x$bblid == ids[i]),]
if (is.element("N",temp[,6]) | is.element("N",temp[,15])) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])}
#try(if (min(temp$mpraxis_rtcr,na.rm=TRUE) < 15) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])})
try(if (max(temp$mpraxis_rtcr,na.rm=TRUE) > 3000) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])})
x[which(x$bblid == ids[i]),] <- temp}  

x <- x[which(x$ntimepoints > 1),]   # best is > 1
x <- x[which(x$timepoint < 5),]     # best is < 6
x <- x[which(is.na(x$cnbAgemonths) == FALSE),]

ids <- unique(x$bblid)
for (i in 1:length(ids)) {
temp <- x[which(x$bblid == ids[i]),]
temp[,17] <- seq(1:nrow(temp))
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

raw_raw <- data.frame(x$mpraxis_rtcr,x$mpraxis_rtcr)
x$mpraxis_rtcr <- scale(x$mpraxis_rtcr)

ACC_r <- scale(winsor(as.matrix(lm(mpraxis_rtcr~Age+Age_Squared+Age_Cubed,data=x)$residuals),trim=0.005))
RT_r <- scale(winsor(as.matrix(lm(mpraxis_rtcr~Age+Age_Squared+Age_Cubed,data=x)$residuals),trim=0.005))

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
if (temp[j,17] == 4) {try(temp[j,21] <- temp[j,4] - temp[which(temp$timepoint == 3),4])}
if (temp[j,17] == 3) {try(temp[j,21] <- temp[j,4] - temp[which(temp$timepoint == 2),4])}
if (temp[j,17] == 2) {try(temp[j,21] <- temp[j,4])}
if (temp[j,17] == 1) {try(temp[j,21] <- NA)}
x[which(x$bblid == ids[i]),] <- temp}}

max_int <- max(x$int,na.rm=TRUE)

int_w_t1 <- x$int
int_w_t1[is.na(int_w_t1) == TRUE] <- max_int
x <- data.frame(x,int_w_t1)

int_sq <- scale(x[,22])^2
int_cub <- scale(x[,22])^3

#regressing out age and practice
ACC_apr <- scale(winsor(as.matrix(lm(mpraxis_rtcr~Age+Age_Squared+Age_Cubed+TP+TP_Squared+TP_Cubed+int_w_t1+int_sq+int_cub,data=x)$residuals),trim=0.005))
RT_apr <- scale(winsor(as.matrix(lm(mpraxis_rtcr~Age+Age_Squared+Age_Cubed+TP+TP_Squared+TP_Cubed+int_w_t1+int_sq+int_cub,data=x)$residuals),trim=0.005))

x$cnbAgemonths <- x$cnbAgemonths/12
x$Time <- x$Time/365.25
Time_Squared <- scale(x$Time)^2
Age <- x$cnbAgemonths
int_ord <- as.factor(quantcut(x$int,3))

x <- data.frame(x[,c(1,3:4,10:11,14:16,16,17,19:22)],ACC_apr,RT_apr,Age,"MPRAXIS",raw_raw)

colnames(x) <- c("	datasetid	","	bblid	","	Time	","	meducation	","	Sex	","	genus	",
"	valid_code	","	ACC_raw	","	RT_raw	","	timepoint	","	ACC_ar	","	RT_ar	","	int	",
"	int_w_t1	","	ACC_apr	","	RT_apr	","	Age	","Test","raw_original_acc","raw_original_RT")

mp <- x


##############################################################################
# TAPPING
##############################################################################

x <- read.csv("CNB Longitudinal/tap_20191113.csv")

x[,4] <- as.numeric(ymd(x[,4])) # date of test
colnames(x)[4] <- "Time"
colnames(x)[3] <- "bblid"
colnames(x)[11] <- "Sex"
x <- x[order(x$bblid, x$timepoint),]

ids <- unique(x$bblid)
for (i in 1:length(ids)) {
temp <- x[which(x$bblid == ids[i]),]
if (is.element("N",temp[,6]) | is.element("N",temp[,15]) | is.element("ctap",temp[,14])) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])}
try(if (min(temp$tap_tot,na.rm=TRUE) < 25) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])})
try(if (max(temp$tap_tot,na.rm=TRUE) > 170) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])})
x[which(x$bblid == ids[i]),] <- temp}  

x <- x[which(x$ntimepoints > 1),]   # best is > 1
x <- x[which(x$timepoint < 5),]     # best is < 6
x <- x[which(is.na(x$cnbAgemonths) == FALSE),]

ids <- unique(x$bblid)
for (i in 1:length(ids)) {
temp <- x[which(x$bblid == ids[i]),]
temp[,21] <- seq(1:nrow(temp))
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

raw_raw <- data.frame(x$tap_tot,x$tap_tot)
x$tap_tot <- scale(x$tap_tot)

ACC_r <- scale(winsor(as.matrix(lm(tap_tot~Age+Age_Squared+Age_Cubed,data=x)$residuals),trim=0.005))
RT_r <- scale(winsor(as.matrix(lm(tap_tot~Age+Age_Squared+Age_Cubed,data=x)$residuals),trim=0.005))

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
if (temp[j,21] == 4) {try(temp[j,25] <- temp[j,4] - temp[which(temp$timepoint == 3),4])}
if (temp[j,21] == 3) {try(temp[j,25] <- temp[j,4] - temp[which(temp$timepoint == 2),4])}
if (temp[j,21] == 2) {try(temp[j,25] <- temp[j,4])}
if (temp[j,21] == 1) {try(temp[j,25] <- NA)}
x[which(x$bblid == ids[i]),] <- temp}}

max_int <- max(x$int,na.rm=TRUE)

int_w_t1 <- x$int
int_w_t1[is.na(int_w_t1) == TRUE] <- max_int
x <- data.frame(x,int_w_t1)

int_sq <- scale(x[,26])^2
int_cub <- scale(x[,26])^3

#regressing out age and practice
ACC_apr <- scale(winsor(as.matrix(lm(tap_tot~Age+Age_Squared+Age_Cubed+TP+TP_Squared+TP_Cubed+int_w_t1+int_sq+int_cub,data=x)$residuals),trim=0.005))
RT_apr <- scale(winsor(as.matrix(lm(tap_tot~Age+Age_Squared+Age_Cubed+TP+TP_Squared+TP_Cubed+int_w_t1+int_sq+int_cub,data=x)$residuals),trim=0.005))

x$cnbAgemonths <- x$cnbAgemonths/12
x$Time <- x$Time/365.25
Time_Squared <- scale(x$Time)^2
Age <- x$cnbAgemonths
int_ord <- as.factor(quantcut(x$int,3))

x <- data.frame(x[,c(1,3:4,10:11,14:15,18,18,21,23:26)],ACC_apr,RT_apr,Age,"TAP",raw_raw)

colnames(x) <- c("	datasetid	","	bblid	","	Time	","	meducation	","	Sex	","	genus	",
"	valid_code	","	ACC_raw	","	RT_raw	","	timepoint	","	ACC_ar	","	RT_ar	","	int	",
"	int_w_t1	","	ACC_apr	","	RT_apr	","	Age	","Test","raw_original_acc","raw_original_RT")

tap <- x


###########################################################################
# COMPILE - this gets written as a "CNB_Longitudinal_Core" file
# write.csv() command omitted below
###########################################################################

x <- rbind(pcet,cpt,nback,cpw,cpf,volt,pvrt,pmat,plOt,er40,medf,adt,mp,tap)
