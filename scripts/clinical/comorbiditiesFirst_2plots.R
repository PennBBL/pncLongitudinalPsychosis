### This script plots percentage of people with each comorbidity AT THE FIRST
### TIME POINT in each of the nine longitudinally defined clinical groups
###
### Ellyn Butler
### December 1, 2020

library(dplyr) # Version 1.0.2
library(ggplot2) # Version 3.3.2
library(ggpubr)

psstat_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/clinical/pnc_longitudinal_diagnosis_n752_202007.csv', stringsAsFactors = TRUE)
diag_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/clinical/n9498_goassess_psych_summary_vars_20131014.csv', stringsAsFactors = TRUE)
demo_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/demographics/baseline/n9498_demo_sex_race_ethnicity_dob.csv', stringsAsFactors = TRUE)

colsforplot <- grep('smry', names(diag_df), value=TRUE)
colsforplot <- colsforplot[!(colsforplot %in% c('smry_prime_pos1', 'smry_prime_tot',
  'smry_prime_pos2', 'smry_psych_overall_rtg', grep('del', names(diag_df), value=TRUE),
  grep('hal', names(diag_df), value=TRUE)))]
colsforplot <- gsub('smry_', '', colsforplot)
colsforplot <- colsforplot[!(colsforplot %in% grep('_cat', colsforplot, value=TRUE))]
colsforplot <- c(colsforplot, 'psy_cat')

diag_df <- diag_df %>%
     mutate_at(grep('smry', names(diag_df), value=TRUE),
      funs(recode(., `0`=0, `1`=0, `2`=0, `3`=0, `4`=1, .default = NaN)))

diag_df <- diag_df[, !(names(diag_df) %in% c('smry_prime_pos1', 'smry_prime_tot',
  'smry_prime_pos2', 'smry_psych_overall_rtg'))]

final_df <- merge(psstat_df, diag_df)
final_df <- merge(final_df, demo_df)

names(final_df) <- gsub('smry_', '', names(final_df))
names(final_df)[names(final_df) == 't1'] <- 'first_diagnosis'
names(final_df)[names(final_df) == 'tfinal2'] <- 'final_diagnosis'
final_df$Female <- recode(final_df$sex, `2`=1, `1`=0)
final_df$sex <- recode(final_df$sex, `2`='Female', `1`='Male')
final_df$White <- recode(final_df$race, `1`=1, .default=0)
final_df$first_diagnosis <- recode(final_df$first_diagnosis, 'other'='OP')
final_df$final_diagnosis <- recode(final_df$final_diagnosis, 'other'='OP')
final_df$t1_tfinal <- recode(final_df$t1_tfinal, 'other_other'='OP_OP',
  'other_TD'='OP_TD', 'other_PS'='OP_PS', 'TD_other'='TD_OP', 'PS_other'='PS_OP')
final_df <- within(final_df, t1_tfinal <- relevel(t1_tfinal, ref='TD_TD'))



### Define percent function
getPercent <- function(i, dataf) {
  first_diag <- as.character(dataf[i, 'first_diagnosis'])
  last_diag <- as.character(dataf[i, 'final_diagnosis'])
  feat <- as.character(dataf[i, 'Feature'])
  (nrow(final_df[final_df$first_diagnosis == first_diag & final_df$final_diagnosis ==
    last_diag & final_df[,feat] == 1  & !is.na(final_df[,feat]),])/nrow(final_df[final_df$first_diagnosis ==
    first_diag & final_df$final_diagnosis == last_diag & !is.na(final_df[,feat]),]))*100
}

### Define N function
getSexN <- function(i, dataf) {
  first_diag <- as.character(dataf[i, 'first_diagnosis'])
  last_diag <- as.character(dataf[i, 'final_diagnosis'])
  feat <- as.character(dataf[i, 'Feature'])
  if (feat == 'Female') { char <- 'F' } else { char <- 'M' }
  paste0(char, ' N = ', nrow(final_df[final_df$first_diagnosis == first_diag & final_df$final_diagnosis ==
    last_diag & final_df$sex == feat & !is.na(final_df$sex),]))
}


##### Get percents
group_df <- read.csv('~/Documents/pncLongitudinalPsychosis/info/firstLastDiags.csv')
group_df <- group_df[, c('First', 'Category')]
group_df <- group_df[!is.na(group_df$First), ]
names(group_df)[names(group_df) == 'First'] <- 'Feature'

final_sum_df <- expand.grid(c('TD', 'OP', 'PS'), c('TD', 'OP', 'PS'),
  colsforplot)
names(final_sum_df) <- c('first_diagnosis', 'final_diagnosis', 'Feature')

final_sum_df$Percent <- sapply(1:nrow(final_sum_df), getPercent, dataf=final_sum_df)

#final_sum_df$first_diagnosis <- recode(final_sum_df$first_diagnosis,
#  'PS'='PS - First Diagnosis', 'OP'='OP - First Diagnosis', 'TD'='TD - First Diagnosis')
#final_sum_df$first_diagnosis <- ordered(final_sum_df$first_diagnosis,
#  c('PS - First Diagnosis', 'OP - First Diagnosis', 'TD - First Diagnosis'))
#final_sum_df$final_diagnosis <- recode(final_sum_df$final_diagnosis,
#  'TD'='TD - Last Diagnosis', 'OP'='OP - Last Diagnosis', 'PS'='PS - Last Diagnosis')
final_sum_df$final_diagnosis <- ordered(final_sum_df$final_diagnosis,
  c('TD', 'OP', 'PS'))

final_sum_df$Feature <- as.character(final_sum_df$Feature)
final_sum_df <- merge(final_sum_df, group_df, by='Feature')


##### Get female and males Ns
ann_text <- expand.grid(c('TD', 'OP', 'PS'), c('TD', 'OP', 'PS'),
  c('Female', 'Male'))
names(ann_text) <- c('first_diagnosis', 'final_diagnosis', 'Feature')

ann_text$lab <- sapply(1:nrow(ann_text), getSexN, dataf=ann_text)

ann_text$first_diagnosis <- recode(ann_text$first_diagnosis,
  'TD'='TD - First Diagnosis', 'OP'='OP - First Diagnosis', 'PS'='PS - First Diagnosis')
ann_text$first_diagnosis <- ordered(ann_text$first_diagnosis,
  c('TD - First Diagnosis', 'OP - First Diagnosis', 'PS - First Diagnosis'))
ann_text$final_diagnosis <- recode(ann_text$final_diagnosis,
  'TD'='TD - Last Diagnosis', 'OP'='OP - Last Diagnosis', 'PS'='PS - Last Diagnosis')
ann_text$final_diagnosis <- ordered(ann_text$final_diagnosis,
  c('TD - Last Diagnosis', 'OP - Last Diagnosis', 'PS - Last Diagnosis'))

sex_labs <- paste0(ann_text[ann_text$Feature == 'Female', 'lab'], ', ',
  ann_text[ann_text$Feature == 'Male', 'lab'])

ann_text <- ann_text[1:9, names(ann_text) != 'Feature']
ann_text$lab <- sex_labs

ann_text$x <- 'dep'
ann_text$y <- 80

final_sum_df$Feature <- gsub('_cat', '', final_sum_df$Feature)
final_sum_df$Feature <- ordered(final_sum_df$Feature, c('ano', 'bul', 'man',
  'add', 'con', 'odd', 'dep', 'gad', 'ocd', 'pan', 'ptd', 'sep', 'soc',
  'agr', 'phb', 'psy')) #Other, Externalizing, Depression, Anxiety, Phobia, Psychosis

final_sum_df$Category <- ordered(final_sum_df$Category, c('Other', 'Externalizing',
  'Depression', 'Anxiety', 'Phobia', 'Psychosis'))

##### Plot
td_baseline_plot <- ggplot(final_sum_df[final_sum_df$first_diagnosis == 'TD',],
    aes(x=Feature, y=Percent, fill=final_diagnosis)) +
  theme_linedraw() + geom_bar(stat='identity', position=position_dodge()) +
  scale_fill_manual(values=c('green3', 'gold', 'red')) +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position='bottom') +
  coord_cartesian(ylim=c(0, 100)) + labs(fill = 'Final Diagnosis') +
  labs(title='TD - First Diagnosis: Comorbidities at First Visit Related to Final Diagnosis') #+
  #geom_text(data=ann_text, mapping=aes(x = x, y = y, label=lab), hjust=-.05, inherit.aes=FALSE)

op_baseline_plot <- ggplot(final_sum_df[final_sum_df$first_diagnosis == 'OP',],
    aes(x=Feature, y=Percent, fill=final_diagnosis)) +
  theme_linedraw() + geom_bar(stat='identity', position=position_dodge()) +
  scale_fill_manual(values=c('green3', 'gold', 'red')) +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position='bottom') +
  coord_cartesian(ylim=c(0, 100)) + labs(fill = 'Final Diagnosis') +
  labs(title='OP - First Diagnosis: Comorbidities at First Visit Related to Final Diagnosis')

ps_baseline_plot <- ggplot(final_sum_df[final_sum_df$first_diagnosis == 'PS',],
    aes(x=Feature, y=Percent, fill=final_diagnosis)) +
  theme_linedraw() + geom_bar(stat='identity', position=position_dodge()) +
  scale_fill_manual(values=c('green3', 'gold', 'red')) +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position='bottom') +
  coord_cartesian(ylim=c(0, 100)) + labs(fill = 'Final Diagnosis') +
  labs(title='PS - First Diagnosis: Comorbidities at First Visit Related to Final Diagnosis')

pdf(file='~/Documents/pncLongitudinalPsychosis/plots/firstComorbid_2plots.pdf', width=14, height=7)
#td_baseline_plot
ggarrange(op_baseline_plot, ps_baseline_plot, ncol=2)
dev.off()
