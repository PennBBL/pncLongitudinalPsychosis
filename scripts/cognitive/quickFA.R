### This script does an exploratory factor analysis of the longitudinal CNB
### data, without accounting for repeated measures. As such, it is important
### that this ultimately be redone.
###
### Ellyn Butler
### August 11, 2020

set.seed(20)

library('dplyr')
library('reshape2')
library('ggplot2')
library('ggpubr')
library('psych')


cnb_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/cognitive/CNB_Longitudinal_Core_11February2020.csv')
cnb_df <- cnb_df[, c('bblid', 'Age', 'timepoint', 'Test', 'ACC_raw')]
cnb_df <- cnb_df[cnb_df$bblid %in% clin_df$bblid,]
cnb_df <- dcast(cnb_df, bblid + Age + timepoint ~ Test, value.var='ACC_raw')

tests <- c('ADT', 'CPF', 'CPT', 'CPW', 'ER40', 'MEDF', 'NBACK', 'PCET',
  'PLOT', 'PMAT', 'PVRT', 'TAP', 'VOLT')
cnb_df[, tests] <- sapply(cnb_df[, tests], scale)
tmp_df <- cnb_df
tmp_df <- tmp_df[!is.na(tmp_df$ADT) & !is.na(tmp_df$CPF) & !is.na(tmp_df$CPT) &
  !is.na(tmp_df$CPW) & !is.na(tmp_df$ER40) & !is.na(tmp_df$MEDF) &
  !is.na(tmp_df$NBACK) & !is.na(tmp_df$PCET) & !is.na(tmp_df$PLOT) &
  !is.na(tmp_df$PMAT) & !is.na(tmp_df$PVRT) & !is.na(tmp_df$TAP) &
  !is.na(tmp_df$VOLT),] # August 11, 2020: EEK. Lose 789 if get rid of rows with any NA

#VSS.scree(tmp_df[,tests])

fanal <- fa(tmp_df[,tests])

tmp_df$MR1 <- factor.scores(tmp_df[, tests], fanal)$scores
tmp_df <- tmp_df[, c('bblid', 'timepoint', 'MR1')]

cnb_df <- merge(cnb_df, tmp_df, all.x=TRUE)

write.csv(cnb_df, '~/Documents/pncLongitudinalPsychosis/data/cognitive/CNB_Longitudinal_Core_Factor_11August2020.csv')
