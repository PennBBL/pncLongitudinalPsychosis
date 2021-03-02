### This script makes table 1 for the clinical paper
###
### Ellyn Butler
### November 30, 2020 - February 11, 2021

library(gtsummary) # Version 1.3.5
library(dbplyr) # Version 1.1.4
library(dplyr) # Version 1.0.2
library(reshape2) # Version 1.4.4

demo_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/demographics/baseline/n9498_demo_sex_race_ethnicity_dob.csv')
#clinical_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/clinical/pnc_longitudinal_diagnosis_n752_202007.csv')
date_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/clinical/pnc_longitudinal_diagnosis_n749_20210112.csv')

# Recalculate ntimepoints (a lot of erroneous zeros)
for (bblid in unique(date_df$bblid)) {
  date_df[date_df$bblid == bblid, 'ntimepoints'] <- length(date_df[date_df$bblid == bblid, 'ntimepoints'])
}

# Create first/last diagnoses df
clinical_df <- date_df[date_df$timepoint == 't1' | date_df$timepoint == paste0('t', date_df$ntimepoints), ]
clinical_df$timepoint <- recode(clinical_df$timepoint, 't1'='t1', 't2'='tfinal2',
  't3'='tfinal2', 't4'='tfinal2', 't5'='tfinal2', 't6'='tfinal2')

clinical_df$diagnosis <- recode(clinical_df$diagnosis, 'psy'='PS')

clinical_df <- reshape2::dcast(clinical_df, bblid ~ timepoint, value.var='diagnosis')
clinical_df$t1_tfinal <- paste(clinical_df$t1, clinical_df$tfinal2, sep='_')

# Merge dataframes
date_df <- merge(date_df, demo_df, by='bblid')
date_df$dob <- as.Date(date_df$dob, '%m/%d/%y')
date_df$dodiagnosis <- as.Date(date_df$dodiagnosis, '%m/%d/%y')
date_df$age <- as.numeric(as.character(((date_df$dodiagnosis - date_df$dob)/365)))


################################### Functions ###################################
# Get age at first and last diagnosis
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

date_df[, c('First Age', 'Final Age')] <- t(sapply(1:nrow(date_df), ageFirstLast))

final_df <- merge(date_df, clinical_df, by='bblid') ##????????????
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

# Positive, Negative, Disorganized, General, GAF, social, role, trauma
sipsgaf_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/clinical/clean_sips_nodup.csv')
names(sipsgaf_df)[names(sipsgaf_df) == 'gaf_c'] <- 'GAF'
sipsgaf_df$AssessmentNumber <- sapply(1:nrow(sipsgaf_df), getAssessmentNumber, dataf=sipsgaf_df)
sipsgaf_df <- sipsgaf_df[sipsgaf_df$AssessmentNumber == 1, ] # Select for first time point

dob_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/demographics/baseline/n9498_demo_sex_race_ethnicity_dob.csv')
corn_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/clinical/Cornblatt_data_all_20210223.csv')
corn_df <- corn_df[corn_df$interview_type %in% c('FP', 'IP', 'YPI'), ]
corn_df <- merge(corn_df, dob_df, by='bblid')
corn_df$dob <- as.Date(corn_df$dob, '%m/%d/%y')
corn_df$date <- as.Date(corn_df$date, '%m/%d/%y')
corn_df$age <- as.numeric(as.character((corn_df$date - corn_df$dob)/365))
corn_df$AssessmentNumber <- sapply(1:nrow(corn_df), getAssessmentNumber, dataf=corn_df)
corn_df <- corn_df[corn_df$AssessmentNumber == 1,] # Select for first time point
names(corn_df)[names(corn_df) == 'cornblatt_gf_role'] <- 'Role'
names(corn_df)[names(corn_df) == 'cornblatt_gf_social'] <- 'Social'
#nrow=820

trauma_df <- read.csv('~/Documents/traumaInformant/data/PNC_GO1_GOASSESSDataArchiveNontext_DATA_2015-07-14_1157.csv')
trauma_df <- trauma_df[trauma_df$interview_type %in% c('AP', 'MP', 'YPI'),]
trauma_df <- trauma_df[,c('proband_bblid', grep('ptd', names(trauma_df), value=TRUE))]
names(trauma_df)[names(trauma_df) == 'proband_bblid'] <- 'bblid'
ptdvars <- c(paste0('ptd00', 1:4), paste0('ptd00', 6:9))
trauma_df[, ptdvars] <- sapply(trauma_df[, ptdvars], na_if, y=9)
trauma_df$Trauma <- rowSums(trauma_df[, ptdvars])


first_df <- merge(final_df, sipsgaf_df, by='bblid', all.x=TRUE)
# Problem bblids: first_df$bblid[duplicated(first_df$bblid)]: 85392  89279  93755 109735 116812
# Two first assessments...



first_df <- merge(first_df, corn_df, by='bblid', all.x=TRUE)
first_df <- merge(first_df, trauma_df, by='bblid', all.x=TRUE)

subsamps <- c('Positive', 'Negative', 'Disorganized', 'General', 'GAF', 'Trauma', 'Social', 'Role')
for (subsamp in subsamps) {
  tmp_df <- first_df[!is.na(first_df[, subsamp]) & !is.na(first_df$Diagnosis), ]
  summary_info <- tmp_df %>%
    dplyr::select('Sex', 'Race', 'Diagnosis', 'First Age') %>%
    tbl_summary(by = Diagnosis)
  gt::gtsave(as_gt(summary_info), paste0('~/Documents/pncLongitudinalPsychosis/results/tableOne_', subsamp, '.html'))
}

summary_info <- final_df %>%
  dplyr::select('Sex', 'Race', 'Diagnosis', 'First Age', 'Final Age') %>%
  tbl_summary(by = Diagnosis)
gt::gtsave(as_gt(summary_info), '~/Documents/pncLongitudinalPsychosis/results/tableOne_full.html')





#
