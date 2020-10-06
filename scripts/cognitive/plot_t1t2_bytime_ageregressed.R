### This script plots male and female trajectories from time 1 to 2 by
### cognitive domain.
###
### Ellyn Butler
### October 5, 2020

set.seed(20)

library('dplyr')
library('ggpubr')
library('reshape2')
library('ggh4x')

cnb_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/cnb_quickFA_impute_sex.csv')

names(cnb_df)[names(cnb_df) == 'sex'] <- 'Sex'

names(cnb_df)[names(cnb_df) == 'EFF_Soln4_MR1'] <- 'SocCog_EFF'
names(cnb_df)[names(cnb_df) == 'EFF_Soln4_MR2'] <- 'Exec_EFF'
names(cnb_df)[names(cnb_df) == 'EFF_Soln4_MR3'] <- 'Mem_EFF'
names(cnb_df)[names(cnb_df) == 'EFF_Soln4_MR4'] <- 'CompCog_EFF'

cnb_df$Sex <- relevel(cnb_df$Sex, 'Male')
#cnb_df$t1_tfinal <- relevel(cnb_df$t1_tfinal, 'TD_TD')

cnb_df$SMSSpeed <- rowMeans(cnb_df[, c('TAP_RT', 'MPRAXIS_RT')])

cogcols <- c('SocCog_EFF', 'Exec_EFF', 'Mem_EFF', 'CompCog_EFF', 'SMSSpeed')

# Regress age out of cognitive scores
cnb_df$agesq <- scale(scale(cnb_df$Age)^2)
cnb_df$agecu <- scale(scale(cnb_df$Age)^3)

cnb_df$SocCog_EFF <- lm(SocCog_EFF ~ Age + agesq + agecu, data=cnb_df)$residuals
cnb_df$Exec_EFF <- lm(Exec_EFF ~ Age + agesq + agecu, data=cnb_df)$residuals
cnb_df$Mem_EFF <- lm(Mem_EFF ~ Age + agesq + agecu, data=cnb_df)$residuals
cnb_df$CompCog_EFF <- lm(CompCog_EFF ~ Age + agesq + agecu, data=cnb_df)$residuals
cnb_df$SMSSpeed <- lm(SMSSpeed ~ Age + agesq + agecu, data=cnb_df)$residuals


# Filter for the first two visits
cnb_df <- cnb_df[cnb_df$Timepoint %in% 1:2, ]
cnb_df <- cnb_df[, c('bblid', 'Age', 'Timepoint', 'Sex', 'first_diagnosis',
  'last_diagnosis', cogcols)]

# Rescale
cnb_df[, cogcols] <- lapply(cnb_df[, cogcols], scale)

# Melt cognitive tests
long_df <- reshape2::melt(cnb_df, c('bblid', 'Age', 'Timepoint', 'Sex',
  'first_diagnosis', 'last_diagnosis'), cogcols)
names(long_df)[names(long_df) == 'variable'] <- 'Test'
names(long_df)[names(long_df) == 'value'] <- 'Score'

long_df$Timepoint <- recode(long_df$Timepoint, `1`='First', `2`='Second')

long_df$first_diagnosis <- recode(long_df$first_diagnosis,
  'PS'='PS - First Diagnosis', 'OP'='OP - First Diagnosis', 'TD'='TD - First Diagnosis')
long_df$first_diagnosis <- ordered(long_df$first_diagnosis,
  c('PS - First Diagnosis', 'OP - First Diagnosis', 'TD - First Diagnosis'))
long_df$last_diagnosis <- recode(long_df$last_diagnosis,
  'TD'='TD - Last Diagnosis', 'OP'='OP - Last Diagnosis', 'PS'='PS - Last Diagnosis')
long_df$last_diagnosis <- ordered(long_df$last_diagnosis,
  c('TD - Last Diagnosis', 'OP - Last Diagnosis', 'PS - Last Diagnosis'))

summary_df <- long_df %>%
        group_by(Timepoint, Sex, Test, first_diagnosis, last_diagnosis) %>%
        summarise(Score = mean(Score))

summary_df$Test <- ordered(summary_df$Test, c('Exec_EFF', 'Mem_EFF',
  'CompCog_EFF', 'SocCog_EFF', 'SMSSpeed'))

#summary_df$Time_Sex_Test <- paste(summary_df$Timepoint, summary_df$Sex, summary_df$Test, sep='_')
summary_df$Time_Test <- paste(summary_df$Timepoint, summary_df$Test, sep='_')
Z <- 1:10
summary_df$Segment <- rep(Z, rep(9, 20))

#https://stackoverflow.com/questions/27005299/ggplot-error-using-linetype-and-group-aesthetics

cnb_plot <- ggplot(summary_df, aes(Test, Score, group=interaction(Sex, Timepoint), colour=Sex, linetype=Timepoint)) +
  theme_linedraw() + facet_nested(first_diagnosis ~ last_diagnosis) + ylim(-.8, .8) +
    	geom_hline(yintercept=0) +
      geom_line(data=summary_df, stat='identity', size=2) +
    	theme(axis.text.x = element_text(angle=90)) +
    	theme(legend.position='bottom', legend.box='vertical', axis.title.x=element_blank(),
    		panel.spacing=unit(.3, 'lines'), axis.text.x = element_text(angle=45, hjust=1)) +
    	scale_color_manual(values=c('blue', 'red'), guide = guide_legend(nrow=1))



pdf(file='~/Documents/pncLongitudinalPsychosis/plots/cnb_t1t2_ageregressed.pdf', width=8, height=7)
cnb_plot
dev.off()
