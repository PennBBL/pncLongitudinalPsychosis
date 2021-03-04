### This script calculates age at various assessments information for the
### participants section
###
### Ellyn Butler
### March 4, 2021

library(dplyr) # Version 1.0.2
library(reshape2) # Version 1.4.4

demo_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/demographics/baseline/n9498_demo_sex_race_ethnicity_dob.csv')
date_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/clinical/pnc_longitudinal_diagnosis_n749_20210112.csv')

# Recalculate ntimepoints (a lot of erroneous zeros)
for (bblid in unique(date_df$bblid)) {
  date_df[date_df$bblid == bblid, 'ntimepoints'] <- length(date_df[date_df$bblid == bblid, 'ntimepoints'])
}

date_df <- merge(date_df, demo_df, by='bblid')
date_df$dob <- as.Date(date_df$dob, '%m/%d/%y')
date_df$dodiagnosis <- as.Date(date_df$dodiagnosis, '%m/%d/%y')
date_df$age <- as.numeric(as.character(((date_df$dodiagnosis - date_df$dob)/365)))

################################### Functions ###################################
ageFirstLast <- function(i) {
  bblid <- date_df[i, 'bblid']
  age_first <- min(date_df[date_df$bblid == bblid, 'age'])
  age_last <- max(date_df[date_df$bblid == bblid, 'age'])
  if (length(age_last) > 0) { # Should always be the case
    c(age_first, age_last)
  } else { c(age_first, NA) }
}

getAssessmentNumber <- function(i, dataf) {
  bblid <- dataf[i, 'bblid']
  ages <- dataf[dataf$bblid == bblid, 'age']
  ages <- sort(unique(ages))
  which(ages == dataf[i, 'age'])
}

#################################################################################

date_df <- date_df[, c('bblid', 'timepoint', 'ntimepoints', 'dodiagnosis',
  'diagnosis', 'sex', 'race', 'ethnicity', 'dob', 'age')]

date_df[, c('first_age', 'final_age')] <- t(sapply(1:nrow(date_df), ageFirstLast))

# Get the age different between first and last assessments
firstlast_df <- date_df[date_df$timepoint == 't1' | date_df$timepoint == paste0('t', date_df$ntimepoints), ]
firstlast_df$timepoint <- recode(firstlast_df$timepoint, 't1'='t1', 't2'='tfinal2',
  't3'='tfinal2', 't4'='tfinal2', 't5'='tfinal2', 't6'='tfinal2')

firstlast_df$diagnosis <- recode(firstlast_df$diagnosis, 'psy'='PS')
firstlast_df <- reshape2::dcast(firstlast_df, bblid + first_age + final_age + ntimepoints~ timepoint, value.var='diagnosis')

mean(firstlast_df$final_age - firstlast_df$first_age)

# Average number of timepoints
mean(firstlast_df$ntimepoints)

# Age range for first and final visits
min(date_df$first_age)
max(date_df$first_age)

min(date_df$final_age)
max(date_df$final_age)















#
