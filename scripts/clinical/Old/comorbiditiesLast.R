### This script plots percentage of people with each comorbidity AT THE FINAL
### TIME POINT in each of the nine longitudinally defined clinical groups
###
### Ellyn Butler
### October 21, 2020 - October 29, 2020

library(dplyr)
library(ggplot2)

psstat_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/clinical/pnc_longitudinal_diagnosis_n752_202007.csv')
diag_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/clinical/pnc_longitudinal_diagnosis_n752_longwdates_202007.csv')
demo_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/demographics/baseline/n9498_demo_sex_race_ethnicity_dob.csv')

diag_df <-  diag_df[diag_df$timepoint == paste0('t', diag_df$ntimepoints), ]
row.names(diag_df) <- 1:nrow(diag_df)

final_df <- merge(psstat_df, diag_df)
final_df <- merge(final_df, demo_df)


names(final_df) <- gsub('dx_', '', names(final_df))
names(final_df)[names(final_df) == 't1'] <- 'first_diagnosis'
names(final_df)[names(final_df) == 'tfinal2'] <- 'last_diagnosis'
final_df$Female <- recode(final_df$sex, `2`=1, `1`=0)
final_df$sex <- recode(final_df$sex, `2`='Female', `1`='Male')
final_df$White <- recode(final_df$race, `1`=1, .default=0)
final_df$first_diagnosis <- recode(final_df$first_diagnosis, 'other'='OP')
final_df$last_diagnosis <- recode(final_df$last_diagnosis, 'other'='OP')
final_df$t1_tfinal <- recode(final_df$t1_tfinal, 'other_other'='OP_OP',
  'other_TD'='OP_TD', 'other_PS'='OP_PS', 'TD_other'='TD_OP', 'PS_other'='PS_OP')
final_df <- within(final_df, t1_tfinal <- relevel(t1_tfinal, ref='TD_TD'))

# Identify bblids PS-TD that have final time point comorbidities (there must be an error here)
final_df[final_df$t1_tfinal == 'PS_TD' & final_df$prodromal_remit == 1, 'bblid'] # 106808
final_df[final_df$t1_tfinal == 'PS_TD' & final_df$sub_abuse == 1, 'bblid'] # 103009
final_df[final_df$t1_tfinal == 'PS_TD' & final_df$sub_abuse_can == 1, 'bblid'] # 103009

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


##### Get percents
final_sum_df <- expand.grid(c('TD', 'OP', 'PS'), c('TD', 'OP', 'PS'),
  names(final_df)[19:40])
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

ann_text$x <- 'moodnos'
ann_text$y <- 80


##### Plot
comorbid_plot <- ggplot(final_sum_df, aes(x=Feature, y=Percent, fill=Feature)) +
  theme_linedraw() + geom_bar(stat='identity') +
  facet_grid(first_diagnosis ~ last_diagnosis) +
  theme(legend.position='none', axis.text.x=element_text(angle=45, hjust=1)) +
  coord_cartesian(ylim=c(0, 100)) +
  #scale_fill_manual(values=c('red', 'pink')) +
  labs(title='Comorbidities at Final Visit') +
  geom_text(data=ann_text, mapping=aes(x = x, y = y, label=lab), hjust=-.05, inherit.aes=FALSE)

pdf(file='~/Documents/pncLongitudinalPsychosis/plots/comorbid3x3.pdf', width=12, height=7)
comorbid_plot
dev.off()
