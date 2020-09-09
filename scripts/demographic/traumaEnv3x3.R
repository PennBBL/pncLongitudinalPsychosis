### This script makes plot of baseline environment and trauma features by
### 3x3 trajectories
###
### Ellyn Butler
### September 9, 2020

library('dplyr')
library('reshape2')
library('ggplot2')
library('ggpubr')
library('fastDummies')
library('broom')
library('stringr')

clin_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/clinical/pnc_longitudinal_diagnosis_n752_202007.csv')
demo_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/demographics/baseline/n1601_demographics_go1_20161212.csv')
trauma_df <- read.csv('~/Documents/traumaInformant/data/PNC_GO1_GOASSESSDataArchiveNontext_DATA_2015-07-14_1157.csv')
env_df <- read.csv('~/Documents/traumaInformant/data/n9498_go1_environment_factor_scores_tymoore_20150909.csv')

trauma_df <- trauma_df[trauma_df$interview_type %in% c('AP', 'MP', 'YPI'),]
trauma_df <- trauma_df[,c('proband_bblid', grep('ptd', names(trauma_df), value=TRUE))]
names(trauma_df)[names(trauma_df) == 'proband_bblid'] <- 'bblid'

final_df <- merge(clin_df, demo_df, by='bblid')
final_df <- merge(final_df, trauma_df, by='bblid')
names(final_df)[names(final_df) == 't1'] <- 'first_diagnosis'
names(final_df)[names(final_df) == 'tfinal2'] <- 'last_diagnosis'
final_df$Female <- recode(final_df$sex, `2`=1, `1`=0)
final_df$White <- recode(final_df$race, `1`=1, .default=0)
final_df$first_diagnosis <- recode(final_df$first_diagnosis, 'other'='OP')
final_df$last_diagnosis <- recode(final_df$last_diagnosis, 'other'='OP')
final_df$t1_tfinal <- recode(final_df$t1_tfinal, 'other_other'='OP_OP',
  'other_TD'='OP_TD', 'other_PS'='OP_PS', 'TD_other'='TD_OP', 'PS_other'='PS_OP')
final_df <- within(final_df, t1_tfinal <- relevel(t1_tfinal, ref='TD_TD'))

ptdvars <- c(paste0('ptd00', 1:4), paste0('ptd00', 6:9))
final_df[, ptdvars] <- sapply(final_df[, ptdvars], na_if, y=9)
final_df <- merge(final_df, env_df)

final_df$sex <- recode(final_df$sex, `2`='Female', `1`='Male')

### Define percent function
getPercent <- function(i, dataf) {
  first_diag <- as.character(dataf[i, 'first_diagnosis'])
  last_diag <- as.character(dataf[i, 'last_diagnosis'])
  feat <- as.character(dataf[i, 'Feature'])
  (nrow(final_df[final_df$first_diagnosis == first_diag & final_df$last_diagnosis ==
    last_diag & final_df[,feat] == 1  & !is.na(final_df[,feat]),])/nrow(final_df[final_df$first_diagnosis ==
    first_diag & final_df$last_diagnosis == last_diag & !is.na(final_df[,feat]),]))*100
}

### Define N function
getSexN <- function(i, dataf) {
  first_diag <- as.character(dataf[i, 'first_diagnosis'])
  last_diag <- as.character(dataf[i, 'last_diagnosis'])
  feat <- as.character(dataf[i, 'Feature'])
  if (feat == 'Female') { char <- 'F' } else { char <- 'M' }
  paste0(char, ' N = ', nrow(final_df[final_df$first_diagnosis == first_diag & final_df$last_diagnosis ==
    last_diag & final_df$sex == feat & !is.na(final_df$sex),]))
}



################################ Sum ################################

final_df$NumTypesTraumas <- rowSums(final_df[, ptdvars])

final_df$Trauma <- recode(final_df$NumTypesTraumas, `0`=0, `1`=1, `2`=1, `3`=1,
  `4`=1, `5`=1, `6`=1, `7`=1)
final_df$BadEnv <- cut(final_df$envHouseholds, quantile(final_df$envHouseholds,
  prob=seq(0, 1, length = 4), type=3))
final_df$BadEnv <- recode(final_df$BadEnv, '(-2.79,-0.431]'=1,
  '(-0.431,0.375]'=0, '(0.375,5.38]'=0)

##### Get percents
final_sum_df <- expand.grid(c('TD', 'OP', 'PS'), c('TD', 'OP', 'PS'),
  c('Trauma', 'BadEnv'))
names(final_sum_df) <- c('first_diagnosis', 'last_diagnosis', 'Feature')

final_sum_df$Percent <- sapply(1:nrow(final_sum_df), getPercent, dataf=final_sum_df)

final_sum_df$first_diagnosis <- recode(final_sum_df$first_diagnosis,
  'PS'='PS - First Diagnosis', 'OP'='OP - First Diagnosis', 'TD'='TD - First Diagnosis')
final_sum_df$first_diagnosis <- ordered(final_sum_df$first_diagnosis,
  c('PS - First Diagnosis', 'OP - First Diagnosis', 'TD - First Diagnosis'))
final_sum_df$last_diagnosis <- recode(final_sum_df$last_diagnosis,
  'TD'='TD - Last Diagnosis', 'OP'='OP - Last Diagnosis', 'PS'='PS - Last Diagnosis')
final_sum_df$last_diagnosis <- ordered(final_sum_df$last_diagnosis,
  c('TD - Last Diagnosis', 'OP - Last Diagnosis', 'PS - Last Diagnosis'))

final_sum_df$Feature <- as.character(final_sum_df$Feature)

##### Get female and males Ns
ann_text <- expand.grid(c('TD', 'OP', 'PS'), c('TD', 'OP', 'PS'),
  c('Female', 'Male'))
names(ann_text) <- c('first_diagnosis', 'last_diagnosis', 'Feature')

ann_text$lab <- sapply(1:nrow(ann_text), getSexN, dataf=ann_text)

ann_text$first_diagnosis <- recode(ann_text$first_diagnosis,
  'TD'='TD - First Diagnosis', 'OP'='OP - First Diagnosis', 'PS'='PS - First Diagnosis')
ann_text$first_diagnosis <- ordered(ann_text$first_diagnosis,
  c('TD - First Diagnosis', 'OP - First Diagnosis', 'PS - First Diagnosis'))
ann_text$last_diagnosis <- recode(ann_text$last_diagnosis,
  'TD'='TD - Last Diagnosis', 'OP'='OP - Last Diagnosis', 'PS'='PS - Last Diagnosis')
ann_text$last_diagnosis <- ordered(ann_text$last_diagnosis,
  c('TD - Last Diagnosis', 'OP - Last Diagnosis', 'PS - Last Diagnosis'))

sex_labs <- paste0(ann_text[ann_text$Feature == 'Female', 'lab'], ', ',
  ann_text[ann_text$Feature == 'Male', 'lab'])

ann_text <- ann_text[1:9, names(ann_text) != 'Feature']
ann_text$lab <- sex_labs

ann_text$x <- 'BadEnv'
ann_text$y <- 80


##### Plot
traumaEnv_plot <- ggplot(final_sum_df, aes(x=Feature, y=Percent, fill=Feature)) +
  theme_linedraw() + geom_bar(stat='identity') +
  facet_grid(first_diagnosis ~ last_diagnosis) +
  theme(legend.position='none') + coord_cartesian(ylim=c(0, 100)) +
  scale_fill_manual(values=c('red', 'pink')) +
  labs(title='Bad Environment or Experience Trauma') +
  geom_text(data=ann_text, mapping=aes(x = x, y = y, label=lab), hjust=-.05, inherit.aes=FALSE)

pdf(file='~/Documents/pncLongitudinalPsychosis/plots/traumaEnv3x3.pdf', width=7, height=7)
traumaEnv_plot
dev.off()
