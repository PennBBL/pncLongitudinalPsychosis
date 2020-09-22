### This script runs GAMMs to see if there are main effects and interactions
### for sex and longitudinal diagnostic labels for the social cognition tests.
### (Sanity check for Kosha)
###
### Ellyn Butler
### September 22, 2020

# September 8, 2020: Bart says look into documentation on 'by' in mgcv to
# understand why factor need to be ordered
#https://stat.ethz.ch/R-manual/R-patched/library/mgcv/html/gam.models.html

set.seed(20)

library('dplyr')
library('reshape2')
library('ggpubr')
library('psych')
library('lme4')
library('gamm4')
library('sjPlot') # Version 2.8.4
library('htmltools')


############################### Social Cog Tests ###############################

clin_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/clinical/pnc_longitudinal_diagnosis_n752_202007.csv')
names(clin_df)[names(clin_df) == 't1'] <- 'first_diagnosis'
names(clin_df)[names(clin_df) == 'tfinal2'] <- 'last_diagnosis'
clin_df$first_diagnosis <- as.character(clin_df$first_diagnosis)
clin_df$last_diagnosis <- as.character(clin_df$last_diagnosis)
clin_df$first_diagnosis <- recode(clin_df$first_diagnosis, 'TD'='TD', 'other'='OP', 'PS'='PS')
clin_df$last_diagnosis <- recode(clin_df$last_diagnosis, 'TD'='TD', 'other'='OP', 'PS'='PS')
clin_df$t1_tfinal <- recode(clin_df$t1_tfinal, 'TD_TD'='TD_TD', 'TD_other'='TD_OP',
  'TD_PS'='TD_PS', 'other_TD'='OP_TD', 'other_other'='OP_OP', 'other_PS'='OP_PS',
  'PS_TD'='PS_TD', 'PS_other'='PS_OP', 'PS_PS'='PS_PS')

demo_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/demographics/baseline/n1601_demographics_go1_20161212.csv')

cnb_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/cognitive/CNB_Longitudinal_Core_11February2020.csv')
cnb_df <- cnb_df[, c('bblid', 'Age', 'timepoint', 'Test', 'ACC_raw', 'RT_raw')]
cnb_df <- cnb_df[cnb_df$bblid %in% clin_df$bblid,]

# Reverse code the RT data (want faster to be higher)
cnb_df$RT_raw <- -cnb_df$RT_raw

cnb_df1 <- dcast(cnb_df, bblid + Age + timepoint ~ Test, value.var='ACC_raw')
names(cnb_df1) <- c('bblid', 'Age', 'Timepoint',
  paste0(names(cnb_df1)[4:length(names(cnb_df1))], '_ACC'))
cnb_df2 <- dcast(cnb_df, bblid + Age + timepoint ~ Test, value.var='RT_raw')
names(cnb_df2) <- c('bblid', 'Age', 'Timepoint',
  paste0(names(cnb_df2)[4:length(names(cnb_df2))], '_RT'))
cnb_df <- merge(cnb_df1, cnb_df2)

cnb_df <- merge(cnb_df, demo_df)
cnb_df$sex <- recode(cnb_df$sex, `2`='Female', `1`='Male')

# Remove TAP_RT, and rename TAP_ACC to TAP_RT
cnb_df <- cnb_df[,!(names(cnb_df) %in% c('TAP_RT', 'MPRAXIS_ACC'))]
names(cnb_df)[names(cnb_df) == 'TAP_ACC'] <- 'TAP_RT'

tests_acc <- c('ADT', 'CPF', 'CPT', 'CPW', 'ER40', 'MEDF', 'NBACK', 'PCET',
  'PLOT', 'PMAT', 'PVRT', 'VOLT')
tests_rt <- c('ADT', 'CPF', 'CPT', 'CPW', 'ER40', 'MEDF', 'NBACK', 'PCET',
  'PLOT', 'PMAT', 'PVRT', 'VOLT', 'TAP', 'MPRAXIS')

cnb_df[, c(paste0(tests_acc, '_ACC'), paste0(tests_rt, '_RT'))] <- sapply(cnb_df[,
  c(paste0(tests_acc, '_ACC'), paste0(tests_rt, '_RT'))], scale)

getDiagnoses <- function(i) {
  bblid <- cnb_df[i, 'bblid']
  c(clin_df[clin_df$bblid == bblid, 'first_diagnosis'], clin_df[clin_df$bblid == bblid, 'last_diagnosis'], as.character(clin_df[clin_df$bblid == bblid, 't1_tfinal']))
}

cnb_df[,c('first_diagnosis', 'last_diagnosis', 't1_tfinal')] <- t(sapply(1:nrow(cnb_df), getDiagnoses))

cnb_df$first_diagnosis <- ordered(cnb_df$first_diagnosis, c('TD', 'OP', 'PS'))
cnb_df$last_diagnosis <- ordered(cnb_df$last_diagnosis, c('TD', 'OP', 'PS'))
cnb_df$t1_tfinal <- relevel(factor(cnb_df$t1_tfinal), ref='TD_TD')




#################################### GAMMs ####################################

cnb_df$sex <- as.factor(as.character(cnb_df$sex))
cnb_df$sex <- relevel(cnb_df$sex, 'Male')
cnb_df$t1_tfinal <- relevel(cnb_df$t1_tfinal, 'TD_TD')

cnb_df$oSex <- ordered(cnb_df$sex, c('Male', 'Female'))

cnb_df$oT1_Tfinal <- ordered(cnb_df$t1_tfinal, c('TD_TD', 'OP_OP', 'OP_PS',
  'OP_TD', 'PS_OP', 'PS_PS', 'PS_TD', 'TD_OP', 'TD_PS'))

tests <- c('ADT_ACC', 'ER40_ACC', 'MEDF_ACC')

test <- 'ADT_ACC'
for (test in tests) {
  test_df <- cnb_df[!is.na(cnb_df[, test]), ]
  row.names(test_df) <- 1:nrow(test_df)

  #### Sex Models
  mod1_sex <- gamm4(as.formula(paste0(test, " ~ sex + s(Age, k=10, bs='cr')")),
    data=test_df, random=~(1|bblid), REML=TRUE)
  mod2_sex <- gamm4(as.formula(paste0(test, " ~ sex + s(Age, k=10, bs='cr') +
    s(Age, by=oSex, k=10, bs='cr')")), data=cnb_df, random=~(1|bblid), REML=TRUE)

  #print(tab_model(mod1_sex$gam, mod2_sex$gam)) Not working
  #save_html(tab_model(mod1_sex$gam, mod2_sex$gam)$page.content, file=paste0('~/Documents/pncLongitudinalPsychosis/results/sex_', test, '.html'))
  #tab_model(mod1_sex$gam, mod2_sex$gam, show.ci=FALSE, show.std="std2", digits.p=5, file=paste0('~/Documents/pncLongitudinalPsychosis/results/sex_', test, '.docx'))
  #tab_model(mod1_sex$gam, mod2_sex$gam, file=paste0('~/Documents/pncLongitudinalPsychosis/results/sex_', test, '.docx'))
  #tab_model(mod1_sex$gam, file=paste0('~/Documents/pncLongitudinalPsychosis/results/sex_', test, '.docx'))
  print(tab_model(mod1_sex$gam, mod2_sex$gam, file=paste0('~/Documents/pncLongitudinalPsychosis/results/sex_', test, '.doc')))

  #### Diagnosis Models
  mod1_diag <- gamm4(as.formula(paste0(test, " ~ t1_tfinal + s(Age, k=10, bs='cr')")),
    data=cnb_df, random=~(1|bblid), REML=TRUE)
  mod2_diag <- gamm4(as.formula(paste0(test, " ~ t1_tfinal + s(Age, k=10, bs='cr') +
    s(Age, by=oT1_Tfinal, k=10, bs='cr')")), data=cnb_df, random=~(1|bblid), REML=TRUE)
  print(tab_model(mod1_diag$gam, mod2_diag$gam, file=paste0('~/Documents/pncLongitudinalPsychosis/results/diag_', test, '.doc')))
}
