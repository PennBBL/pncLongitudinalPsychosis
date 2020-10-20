#clear workspace
setwd()
#start code
library(ggplot2)
#install.packages("devtools")
library(devtools)
#install_github("easyGgplot2", "kassambara")
library(easyGgplot2)
library(pROC)
library(doBy)
library(MASS)
library(reshape2)
library(mgcv)
library(gamm4)
library(RLRsim)
library(car)
#install.packages("psych")
library(psych)
#install.packages("car")
library(car)
#install.packages("visreg")
library(visreg)


diagnosis <-read.csv("~/Documents/pncLongitudinalPsychosis/data/clinical/pnc_longitudinal_diagnosis_n752_202007.csv")
time2dates <-read.csv("~/Documents/pncLongitudinalPsychosis/data/clinical/pnc_longitudinal_diagnosis_n752_longwdates_202007.csv")
Diagnosiswgroups1 <- merge(diagnosis, time2dates, by.x ="bblid", all.x = FALSE)


Diagnosiswgroups2 <- Diagnosiswgroups1[,c(1,3,10,11,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40)]
View(Diagnosiswgroups2)
Diagnosiswgroups2 %>% distinct("bblid", .keep_all = TRUE)


library(dplyr)
Diagnosiswgroups<-distinct(Diagnosiswgroups2,bblid,dx_moodnos, .keep_all = TRUE)
Diagnosiswgroups<-distinct(Diagnosiswgroups2,bblid,dx_mdd, .keep_all = TRUE)
Diagnosiswgroups<-distinct(Diagnosiswgroups2,bblid,dx_bp1, .keep_all = TRUE)
Diagnosiswgroups<-distinct(Diagnosiswgroups2,bblid,dx_bpoth, .keep_all = TRUE)
Diagnosiswgroups<-distinct(Diagnosiswgroups2,bblid,dx_adhd, .keep_all = TRUE)
Diagnosiswgroups<-distinct(Diagnosiswgroups2,bblid,dx_anx, .keep_all = TRUE)
Diagnosiswgroups<-distinct(Diagnosiswgroups2,bblid,dx_ptsd, .keep_all = TRUE)
Diagnosiswgroups<-distinct(Diagnosiswgroups2,bblid,dx_cogdis, .keep_all = TRUE)
Diagnosiswgroups<-distinct(Diagnosiswgroups2,bblid,dx_other, .keep_all = TRUE)
Diagnosiswgroups<-distinct(Diagnosiswgroups2,bblid,dx_BrderPD, .keep_all = TRUE)
Diagnosiswgroups<-distinct(Diagnosiswgroups2,bblid,dx_sub_dep_can, .keep_all = TRUE)
Diagnosiswgroups<-distinct(Diagnosiswgroups2,bblid,dx_sub_dep_alc, .keep_all = TRUE)
Diagnosiswgroups<-distinct(Diagnosiswgroups2,bblid,dx_sub_dep_oth, .keep_all = TRUE)
Diagnosiswgroups<-distinct(Diagnosiswgroups2,bblid,dx_sub_abuse_can, .keep_all = TRUE)
Diagnosiswgroups<-distinct(Diagnosiswgroups2,bblid,dx_sub_abuse_alc, .keep_all = TRUE)
Diagnosiswgroups<-distinct(Diagnosiswgroups2,bblid,dx_sub_abuse_oth, .keep_all = TRUE)
Diagnosiswgroups<-distinct(Diagnosiswgroups2,bblid,dx_sub_dep, .keep_all = TRUE)
Diagnosiswgroups<-distinct(Diagnosiswgroups2,bblid,dx_sub_abuse, .keep_all = TRUE)



DiagnosiswgroupsTD_TD<-Diagnosiswgroups[ which(Diagnosiswgroups$t1_tfinal=='TD_TD'),]

Group <- c('dx_moodnos', 'dx_mdd','dx_bp1','dx_bpoth','dx_adhd','dx_anx','dx_ptsd','dx_cogdis','dx_other','dx_BrderPD','dx_sub_dep_can','dx_sub_dep_alc','dx_sub_dep_oth','dx_sub_abuse_can','dx_sub_abuse_alc','dx_sub_abuse_oth','dx_sub_dep','dx_sub_abuse')
Number <- c((sum(DiagnosiswgroupsTD_TD$dx_moodnos, na.rm=T)),sum(DiagnosiswgroupsTD_TD$dx_mdd, na.rm=T),sum(DiagnosiswgroupsTD_TD$dx_bp1, na.rm=T),sum(DiagnosiswgroupsTD_TD$dx_bpoth, na.rm=T),
            sum(DiagnosiswgroupsTD_TD$dx_adhd, na.rm=T),sum(DiagnosiswgroupsTD_TD$dx_anx, na.rm=T),
            sum(DiagnosiswgroupsTD_TD$dx_ptsd, na.rm=T),sum(DiagnosiswgroupsTD_TD$dx_cogdis, na.rm=T),sum(DiagnosiswgroupsTD_TD$dx_other, na.rm=T),
            sum(DiagnosiswgroupsTD_TD$dx_BrderPD, na.rm=T),sum(DiagnosiswgroupsTD_TD$dx_sub_dep_can, na.rm=T),
            sum(DiagnosiswgroupsTD_TD$dx_sub_dep_alc, na.rm=T),sum(DiagnosiswgroupsTD_TD$dx_sub_dep_oth, na.rm=T),
            sum(DiagnosiswgroupsTD_TD$dx_sub_abuse_can, na.rm=T),sum(DiagnosiswgroupsTD_TD$dx_sub_abuse_alc, na.rm=T),
            sum(DiagnosiswgroupsTD_TD$dx_sub_abuse_oth, na.rm=T),sum(DiagnosiswgroupsTD_TD$dx_sub_dep, na.rm=T),
            sum(DiagnosiswgroupsTD_TD$dx_sub_abuse, na.rm=T))
Total <- c(length(unique(DiagnosiswgroupsTD_TD[["bblid"]])))
Percentage <- c((Number/Total)*100)
Diagnosis1 <- c("TD")
Diagnosis2 <- c("TD")
TD_TD<- data.frame(Diagnosis1,Diagnosis2,Group,Number,Total,Percentage)
View(TD_TD)

library(concatenate)
newframe <- rbind(OP_PS,OP_OP,OP_TD,PS_PS,PS_TD,PS_OP,TD_OP,TD_PS,TD_TD)
View(newframe)

###makeplot

library(ggplot2)
PS_TDplot <- ggplot(PS_TD, aes(x=reorder(Group, Percentage), y=Percentage)) +
  geom_bar(stat = "identity", colour = "black", width = .5) + 
  scale_x_discrete("") +
  scale_y_discrete("Percentage of PS_TD", limits = factor(0:50)) +
  theme_linedraw() +
  theme(panel.grid.major.y = element_blank())

PS_TDplot

library(ggplot2)
newframeplot <- ggplot(newframe, aes(x=reorder(Group,Percentage), y=Percentage)) +
  geom_bar(stat = "identity", colour = "black", width = .5) + 
  scale_x_discrete("", labels = c("Cogdis",'BrderPD','SubAbuse_Other','BP1','Other','PTSD','SubDep_Other','BP_Other','SubDep_Alc','SubAbuse_Can','SubAbuse_Alc','Anxiety','Moodnos','SubDep_Can',
                   'SubAbuse','SubDep','ADHD','MDD')) +
  scale_y_discrete("Percentage", limits = factor(0:50), breaks=c(1,5,10,15,20,25,30,35,40,45,50)) +
  theme_linedraw() +
  theme(panel.grid.major.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 45)) 
  
####

newframeplot + facet_grid(Diagnosis1 ~ Diagnosis2)
pdf("new")




###PS_TD

DiagnosiswgroupsTD_TD<-Diagnosiswgroups[ which(Diagnosiswgroups$t1_tfinal=='PS_PS'),]
Group <- c('dx_moodnos', 'dx_mdd','dx_bp1','dx_bpoth','dx_adhd','dx_anx','dx_ptsd','dx_cogdis','dx_other','dx_BrderPD','dx_sub_dep_can','dx_sub_dep_alc','dx_sub_dep_oth','dx_sub_abuse_can','dx_sub_abuse_alc','dx_sub_abuse_oth','dx_sub_dep','dx_sub_abuse')
Number <- c((sum(DiagnosiswgroupsTD_TD$dx_moodnos, na.rm=T)),sum(DiagnosiswgroupsTD_TD$dx_mdd, na.rm=T),sum(DiagnosiswgroupsTD_TD$dx_bp1, na.rm=T),sum(DiagnosiswgroupsTD_TD$dx_bpoth, na.rm=T),
            sum(DiagnosiswgroupsTD_TD$dx_adhd, na.rm=T),sum(DiagnosiswgroupsTD_TD$dx_anx, na.rm=T),
            sum(DiagnosiswgroupsTD_TD$dx_ptsd, na.rm=T),sum(DiagnosiswgroupsTD_TD$dx_cogdis, na.rm=T),sum(DiagnosiswgroupsTD_TD$dx_other, na.rm=T),
            sum(DiagnosiswgroupsTD_TD$dx_BrderPD, na.rm=T),sum(DiagnosiswgroupsTD_TD$dx_sub_dep_can, na.rm=T),
            sum(DiagnosiswgroupsTD_TD$dx_sub_dep_alc, na.rm=T),sum(DiagnosiswgroupsTD_TD$dx_sub_dep_oth, na.rm=T),
            sum(DiagnosiswgroupsTD_TD$dx_sub_abuse_can, na.rm=T),sum(DiagnosiswgroupsTD_TD$dx_sub_abuse_alc, na.rm=T),
            sum(DiagnosiswgroupsTD_TD$dx_sub_abuse_oth, na.rm=T),sum(DiagnosiswgroupsTD_TD$dx_sub_dep, na.rm=T),
            sum(DiagnosiswgroupsTD_TD$dx_sub_abuse, na.rm=T))
Total <- c(length(unique(DiagnosiswgroupsTD_TD[["bblid"]])))
Percentage <- c((Number/Total)*100)
PS_PS<- data.frame(Group,Number,Total,Percentage)
View(PS_PS)









##PS_PS plot done 
##PS_TD plot done
##PS_other 




###normal orientation
p <- ggplot(testPS_TD, aes(x=reorder(Group, -Percentage), y=Percentage, drop(y=0))) +
      geom_bar(stat = "identity", colour = "black", width = .5) +
  scale_x_discrete("TD_PS")

####flip orientation w/ same order 
p <- ggplot(testPS_TD, aes(x=Group, y=Percentage, drop(y=0))) +
  geom_bar(stat = "identity", colour = "black", width = .5) +
  scale_x_discrete("Comorbidities") +
  scale_y_discrete("Percentage of TD_PS", limits = factor(0:10)) +
  coord_flip() 

#### flip orientation w/ highest to lowest
TD_PSplot <- ggplot(testPS_TD, aes(x=reorder(Group, Percentage), y=Percentage, drop(y=0))) +
  geom_bar(stat = "identity", colour = "black", width = .5) +
  scale_x_discrete("Comorbidities") + 
  scale_y_discrete("Percentage of TD_PS", limits = factor(0:10)) +
  coord_flip()
TD_PSplot


### same code for TD_other
DiagnosiswgroupsTD_TD<-Diagnosiswgroups[ which(Diagnosiswgroups$t1_tfinal=='PS_PS'),]
View(DiagnosiswgroupsTD_TD)

###get numbers for each group
Group <- c('dx_moodnos', 'dx_mdd','dx_bp1','dx_bpoth','dx_adhd','dx_anx','dx_ptsd','dx_cogdis','dx_other','dx_BrderPD','dx_sub_dep_can','dx_sub_dep_alc','dx_sub_dep_oth','dx_sub_abuse_can','dx_sub_abuse_alc','dx_sub_abuse_oth','dx_sub_dep','dx_sub_abuse')
Number <- c((sum(DiagnosiswgroupsTD_TD$dx_moodnos, na.rm=T)),sum(DiagnosiswgroupsTD_TD$dx_mdd, na.rm=T),sum(DiagnosiswgroupsTD_TD$dx_bp1, na.rm=T),sum(DiagnosiswgroupsTD_TD$dx_bpoth, na.rm=T),
            sum(DiagnosiswgroupsTD_TD$dx_adhd, na.rm=T),sum(DiagnosiswgroupsTD_TD$dx_anx, na.rm=T),
            sum(DiagnosiswgroupsTD_TD$dx_ptsd, na.rm=T),sum(DiagnosiswgroupsTD_TD$dx_cogdis, na.rm=T),sum(DiagnosiswgroupsTD_TD$dx_other, na.rm=T),
            sum(DiagnosiswgroupsTD_TD$dx_BrderPD, na.rm=T),sum(DiagnosiswgroupsTD_TD$dx_sub_dep_can, na.rm=T),
            sum(DiagnosiswgroupsTD_TD$dx_sub_dep_alc, na.rm=T),sum(DiagnosiswgroupsTD_TD$dx_sub_dep_oth, na.rm=T),
            sum(DiagnosiswgroupsTD_TD$dx_sub_abuse_can, na.rm=T),sum(DiagnosiswgroupsTD_TD$dx_sub_abuse_alc, na.rm=T),
            sum(DiagnosiswgroupsTD_TD$dx_sub_abuse_oth, na.rm=T),sum(DiagnosiswgroupsTD_TD$dx_sub_dep, na.rm=T),
            sum(DiagnosiswgroupsTD_TD$dx_sub_abuse, na.rm=T))
Total <- c(length(unique(DiagnosiswgroupsTD_TD[["bblid"]])))
Percentage <- c((Number/Total)*100)
PS_PS <- data.frame(Group,Number,Total,Percentage)
View(TD_other)





graph1 <- grid.arrange(PS_TDplot, PS_PSplot, PS_PSplot, nrow=1)
graph2 <- grid.arrange(PS_PSplot, PS_TDplot, PS_TDplot, nrow=1)
graph3 <- grid.arrange(PS_TDplot, PS_PSplot, PS_TDplot, nrow=1)
testfacet <- grid.arrange(graph1, graph2, graph3, nrow=2)
                       
                       grid.arrange(TD_OPplot, TD_PSplot, nrow=1)
testfacet



graph1





Diagnosiswgroups<-Diagnosiswgroups[ which(Diagnosiswgroups$t1_tfinal != 'TD_other'),]
Diagnosiswgroups<-Diagnosiswgroups[ which(Diagnosiswgroups$t1_tfinal != 'TD_TD'),]
Diagnosiswgroups<-Diagnosiswgroups[ which(Diagnosiswgroups$t1_tfinal != 'PS_other'),]
Diagnosiswgroups<-Diagnosiswgroups[ which(Diagnosiswgroups$t1_tfinal != 'other_other'),]
Diagnosiswgroups<-Diagnosiswgroups[ which(Diagnosiswgroups$t1_tfinal != 'other_TD'),]
Diagnosiswgroups<-Diagnosiswgroups[ which(Diagnosiswgroups$t1_tfinal != 'PS_PS'),]
Diagnosiswgroups<-Diagnosiswgroups[ which(Diagnosiswgroups$t1_tfinal != 'PS_other'),]
Diagnosiswgroups<-Diagnosiswgroups[ which(Diagnosiswgroups$t1_tfinal != 'PS_TD'),]
Diagnosiswgroups<-Diagnosiswgroups[ which(Diagnosiswgroups$t1_tfinal != 'other_PS'),]



  
