### This script makes table 1 for the clinical paper
###
### Ellyn Butler
### November 30, 2020

library('dbplyr') # Version 1.1.4
library('gtsummary') # Version 1.3.5
library('dplyr') # Version 1.0.2

demo_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/demographics/baseline/n9498_demo_sex_race_ethnicity_dob.csv')
clinical_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/clinical/pnc_longitudinal_diagnosis_n752_202007.csv')
date_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/clinical/pnc_longitudinal_diagnosis_n752_longwdates_202007.csv')

date_df <- merge(date_df, demo_df, by='bblid')
date_df$dob <- as.Date(date_df$dob, '%m/%d/%y')
date_df$DODIAGNOSIS <- as.Date(date_df$DODIAGNOSIS, '%m/%d/%y')
date_df$age <- (date_df$DODIAGNOSIS - date_df$dob)/365


# Get age at first and last diagnosis
ageFirstLast <- function(i) {
  bblid <- date_df[i, 'bblid']
  age_first <- min(date_df[date_df$bblid == bblid, 'age'])
  age_last <- max(date_df[date_df$bblid == bblid, 'age'])
  if (length(age_last) > 0) { # Should always be the case
    c(age_first, age_last)
  } else { c(age_first, NA) }
}

date_df[, c('First Age', 'Final Age')] <- t(sapply(1:nrow(date_df), ageFirstLast))

final_df <- merge(date_df, clinical_df, by='bblid')
final_df <- final_df[final_df$timepoint == 't1', ]
row.names(final_df) <- 1:nrow(final_df)

# Recode variables
final_df$Sex <- recode(final_df$sex, `1`='Male', `2`='Female')
final_df$Race <- recode(final_df$race, `1`='Caucasian', `2`='African American',
  `3`='US India/Alaska Native', `4`='Asian', `5`='More Than One Race')
final_df$Diagnosis <- recode(final_df$t1_tfinal, 'other_TD'='OP-TD',
  'other_other'='OP-OP', 'other_PS'='OP-PS', 'TD_other'='TD-OP', 'PS_other'='PS-OP',
  'TD_TD'='TD-TD', 'TD_PS'='TD-PS', 'PS_TD'='PS-TD', 'PS_PS'='PS-PS')
final_df$Diagnosis <- ordered(final_df$Diagnosis, c('TD-TD', 'TD-OP', 'TD-PS',
  'OP-TD', 'OP-OP', 'OP-PS', 'PS-TD', 'PS-OP', 'PS-PS'))

summary_info <- final_df %>%
  dplyr::select('Sex', 'Race', 'Diagnosis', 'First Age', 'Final Age') %>%
  tbl_summary(by = Diagnosis)
