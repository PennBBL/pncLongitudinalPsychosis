### Test if TD-PS is higher on psychosis factor score at baseline than TD-TD
### and TD-OP.
###
### Ellyn Butler
### October 5, 2020

library('dplyr')
library('reshape2')

clin_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/clinical/pnc_longitudinal_diagnosis_n752_202007.csv')
demo_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/demographics/baseline/n1601_demographics_go1_20161212.csv')
psychosis_df <- read.csv('~/Documents/predLongLabels/data/n9498_goassess_itemwise_corrtraits_scores_20161219.csv')

final_df <- merge(clin_df, demo_df, by='bblid')
final_df <- merge(final_df, psychosis_df, by='bblid')
names(final_df)[names(final_df) == 't1'] <- 'first_diagnosis'
names(final_df)[names(final_df) == 'tfinal2'] <- 'last_diagnosis'
final_df$Female <- recode(final_df$sex, `2`=1, `1`=0)
final_df$White <- recode(final_df$race, `1`=1, .default=0)
final_df$first_diagnosis <- recode(final_df$first_diagnosis, 'other'='OP')
final_df$last_diagnosis <- recode(final_df$last_diagnosis, 'other'='OP')
final_df$t1_tfinal <- recode(final_df$t1_tfinal, 'other_other'='OP_OP',
  'other_TD'='OP_TD', 'other_PS'='OP_PS', 'TD_other'='TD_OP', 'PS_other'='PS_OP')
final_df <- within(final_df, t1_tfinal <- relevel(t1_tfinal, ref='TD_TD'))


t.test(final_df[final_df$t1_tfinal == 'TD_PS', 'psychosis_corrtraitsv2'],
  final_df[final_df$t1_tfinal == 'TD_TD', 'psychosis_corrtraitsv2'])

t.test(final_df[final_df$t1_tfinal == 'TD_PS', 'psychosis_corrtraitsv2'],
  final_df[final_df$t1_tfinal == 'TD_OP', 'psychosis_corrtraitsv2'])
