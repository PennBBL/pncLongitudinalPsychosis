### This script plots age at worst experience by longitudinal clinical labels
###
### Ellyn Butler
### August 12, 2020

library('dplyr')
library('reshape2')
library('ggplot2')
library('ggpubr')

clin_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/clinical/pnc_longitudinal_diagnosis_n752_202007.csv')

trauma_df <- read.csv('~/Documents/traumaInformant/data/PNC_GO1_GOASSESSDataArchiveNontext_DATA_2015-07-14_1157.csv')
trauma_df <- trauma_df[trauma_df$interview_type %in% c('AP', 'MP', 'YPI'),]
trauma_df <- trauma_df[,c('proband_bblid', grep('ptd', names(trauma_df), value=TRUE))]
names(trauma_df)[names(trauma_df) == 'proband_bblid'] <- 'bblid'

final_df <- merge(clin_df, trauma_df, by='bblid')
final_df$NumTypesTraumas <- rowSums(final_df[, c(paste0('ptd00', 1:4), paste0('ptd00', 6:9))])
final_df <- final_df[final_df$NumTypesTraumas > 0,]
