setwd("~/Documents/JENNA pnc longituidinal psychosis ")
library(dbplyr)
library(gtsummary)
install.packages("dplyr")   
demo <-read.csv("n9498_demo_sex_race_ethnicity_dob.csv")
diagnosis <-read.csv("pnc_longitudinal_diagnosis_n752_202007.csv")
time2dates <-read.csv("pnc_longitudinal_diagnosis_n752_longwdates_202007.csv")
demo2 <-read.csv("n1601_demographics_go1_20161212.csv")
finalages <-read.csv("ageatassessments.csv")

maindata1 <-merge(demo, diagnosis, by="bblid")
totaldata <-merge(maindata1, finalages, by="bblid")

###time2 <- merge(totaldata, time2dates, by="bblid")
####write.csv(time2, "time2DOD.csv")

###totaltime2 <- read.csv("time2DOD.csv")
###View(totaltime2)

### get rid of columns not needed
##importantinfo <- maindatafinalage[,c(1,2,3,4,5,7,14,15,16)]
###write.csv(importantinfo, "importantinfo.csv")
###View(importantinfo)

names(totaldata)[2] <- "Sex"
names(totaldata)[3] <- "Race"
names(totaldata)[4] <- "Ethinicity"
totaldata$Race[totaldata$Race=="1"] <- "Caucasian"
totaldata$Race[totaldata$Race=="2"] <- "African American"
totaldata$Race[totaldata$Race=="3"] <- "US India/Alaska Native"
totaldata$Race[totaldata$Race=="4"] <- "Asian"
totaldata$Race[totaldata$Race=="5"]<- "More Than One Race"
totaldata$Sex[totaldata$Sex=="1"]<- "Male"
totaldata$Sex[totaldata$Sex=="2"]<- "Female"
totaldata$t1[totaldata$t1=="other"]<-"OP"
totaldata$t1[totaldata$t1_tfinal=="other"]<-"OP"
totaldata$t1_tfinal[totaldata$t1_tfinal=="TD_other"]<-"TD_OP"
totaldata$t1_tfinal[totaldata$t1_tfinal=="PS_other"]<-"PS_OP"
totaldata$t1_tfinal[totaldata$t1_tfinal=="other_other"]<-"OP_OP"
totaldata$t1_tfinal[totaldata$t1_tfinal=="other_PS"]<-"OP_PS"
totaldata$t1_tfinal[totaldata$t1_tfinal=="other_TD"]<-"OP_TD"
View(totaldata)

testtest4 <- totaldata %>% 
  dplyr::select(Sex, Race, t1_tfinal, ageAtTime1, ageAtFinalAssess)
testtest4 %>% tbl_summary(by = t1_tfinal) 

######unused code below 

names(importantinfo)[2] <- "Sex"
names(importantinfo)[3] <- "Race"
names(importantinfo)[4] <- "Ethinicity"
importantinfo$Race[importantinfo$Race=="1"]<-"Cauc"
importantinfo$Race[importantinfo$Race=="2"]<-"AA"
importantinfo$Race[importantinfo$Race=="3"]<-"Other"
importantinfo$Race[importantinfo$Race=="4"]<-"Other"
importantinfo$Race[importantinfo$Race=="5"]<-"Other"
importantinfo$Sex[importantinfo$Sex=="1"]<-"Male"
importantinfo$Sex[importantinfo$Sex=="2"]<-"Female"
View(importantinfo)

time2 <- merge(demo, time2dates, by="bblid")
dateofdiag <- time2[,c(1,5,6)]
View(dateofdiag)
age_calc(dob, enddate=DO)

####testtest5 <- totaldata %>% 
  dplyr::select(Sex, Race, t1, atAtTime1)
testtest5 %>% 
  tbl_summary(by = t1) %>% add_p()


####testtest6 <- importantinfo %>% 
  dplyr::select(Sex, Race, t1_tfinal, AgeatAsessment)
testtest6 %>% 
  tbl_summary(by = t1_tfinal) 

####teststuff
teststuff <- merge(totaldata,totaltime2, by="bblid", no.dups = FALSE)






