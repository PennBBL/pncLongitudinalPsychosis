### This script tests if the time difference between first and last assessments
### explains any variance in final PS status.
###
### Ellyn Butler
### October 5, 2020 - October 19, 2020

library('dplyr')
library('sjPlot')
library('reshape2')
library('fastDummies')
library('lubridate')

clin_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/clinical/pnc_longitudinal_diagnosis_n752_longwdates_202007.csv')
demo_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/demographics/baseline/n9498_demo_sex_race_ethnicity_dob.csv')
env_df <- read.csv('~/Documents/traumaInformant/data/n9498_go1_environment_factor_scores_tymoore_20150909.csv')
trauma_df <- read.csv('~/Documents/traumaInformant/data/PNC_GO1_GOASSESSDataArchiveNontext_DATA_2015-07-14_1157.csv')
diag_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/clinical/pnc_longitudinal_diagnosis_n752_202007.csv')

final_df <- merge(clin_df, demo_df, by='bblid')
final_df <- merge(final_df, diag_df, by='bblid')

names(final_df)[names(final_df) == 't1'] <- 'first_diagnosis'
names(final_df)[names(final_df) == 'tfinal2'] <- 'last_diagnosis'
final_df$Female <- recode(final_df$sex, `2`=1, `1`=0)
final_df$White <- recode(final_df$race, `1`=1, .default=0)
final_df$first_diagnosis <- recode(final_df$first_diagnosis, 'other'='OP')
final_df$last_diagnosis <- recode(final_df$last_diagnosis, 'other'='OP')
final_df$t1_tfinal <- recode(final_df$t1_tfinal, 'other_other'='OP_OP',
  'other_TD'='OP_TD', 'other_PS'='OP_PS', 'TD_other'='TD_OP', 'PS_other'='PS_OP')
final_df <- within(final_df, t1_tfinal <- relevel(t1_tfinal, ref='TD_TD'))
final_df <- within(final_df, first_diagnosis <- relevel(first_diagnosis, ref='TD'))

final_df <- merge(final_df, env_df)
final_df$PS_final <- recode(final_df$tfinal, 'other'='No', 'PS'='Yes',
  'psy'='Yes', 'TD'='No')

# Calculate ages
final_df$DODIAGNOSIS <- as.Date(final_df$DODIAGNOSIS, '%m/%d/%y')
final_df$dob <- as.Date(final_df$dob, '%m/%d/%y')
final_df$Age <- lubridate::time_length(difftime(final_df$DODIAGNOSIS, final_df$dob), 'years')


getFirstLastAge <- function(i) {
  bblid <- final_df[i, 'bblid']
  age_first <- min(final_df[final_df$bblid == bblid, 'Age'])
  age_last <- max(final_df[final_df$bblid == bblid, 'Age'])
  c(age_first, age_last)
}

final_df[, c('age_first', 'age_last')] <- t(sapply(1:nrow(final_df), getFirstLastAge))

final_df$timediff_t1_tfinal <- final_df$age_last - final_df$age_first

#final_df <- merge(final_df, clin_df)

trauma_df <- trauma_df[trauma_df$interview_type %in% c('AP', 'MP', 'YPI'),]
trauma_df <- trauma_df[,c('proband_bblid', grep('ptd', names(trauma_df), value=TRUE))]
names(trauma_df)[names(trauma_df) == 'proband_bblid'] <- 'bblid'

final_df <- merge(final_df, trauma_df)

ptdvars <- c(paste0('ptd00', 1:4), paste0('ptd00', 6:9))
final_df[, ptdvars] <- sapply(final_df[, ptdvars], na_if, y=9)
final_df$num_type_trauma <- rowSums(final_df[, ptdvars])


##### Filter for only the first time point
final_df <- final_df[final_df$timepoint == 't1',]

######### Does # of trauma types predict final time point PS status? #########

mod0 <- glm(PS_final ~ Age, family='binomial', data=final_df)

mod1 <- glm(PS_final ~ Age + timediff_t1_tfinal, family='binomial', data=final_df)

mod2 <- glm(PS_final ~ Age + timediff_t1_tfinal + first_diagnosis, family='binomial', data=final_df)

mod3 <- glm(PS_final ~ Age + timediff_t1_tfinal + first_diagnosis + num_type_trauma, family='binomial', data=final_df)

#Answer: Yes

# Does baseline diagnosis moderate the relationship between final time
# point PS status and number of types of traumas such that there is a stronger
# relationship between number of trauma types and final PS status among those
# who were PS at baseline than those who were OP or TD?

mod4 <- glm(PS_final ~ Age + timediff_t1_tfinal + num_type_trauma*first_diagnosis, family='binomial', data=final_df)

#Answer: No

all_models_trauma <- tab_model(mod0, mod1, mod2, mod3, mod4)
