### This script runs GAMMs to see if there are main effects and interactions
### for sex and longitudinal SIPS scores
### score.
###
### Ellyn Butler
### October 6, 2020

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
library('sjPlot')

clin_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/clinical/clean_sips.csv')
diag_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/clinical/pnc_longitudinal_diagnosis_n752_202007.csv')
diag_df$t1_tfinal <- recode(diag_df$t1_tfinal, 'TD_TD'='TD_TD', 'TD_other'='TD_OP',
  'TD_PS'='TD_PS', 'other_TD'='OP_TD', 'other_other'='OP_OP', 'other_PS'='OP_PS',
  'PS_TD'='PS_TD', 'PS_other'='PS_OP', 'PS_PS'='PS_PS')

clin_df <- merge(clin_df, diag_df, by='bblid')

clin_df$sex <- relevel(clin_df$sex, 'Male')
clin_df$t1_tfinal <- relevel(clin_df$t1_tfinal, 'TD_TD')

clin_df$oT1_Tfinal <- ordered(clin_df$t1_tfinal, c('TD_TD', 'OP_OP', 'OP_PS',
  'OP_TD', 'PS_OP', 'PS_PS', 'PS_TD', 'TD_OP', 'TD_PS'))

names(clin_df)[names(clin_df) == 'gaf_c'] <- 'GAF_Current'
names(clin_df)[names(clin_df) == 'gaf_h'] <- 'GAF_Highest'
names(clin_df)[names(clin_df) == 'age'] <- 'Age'


################################## GAMMs ##################################

############## Positive ##############
mod1_diag <- gamm4(Positive ~ t1_tfinal + s(Age, k=4, bs='cr'), data=clin_df,
  random=~(1|bblid), REML=TRUE)
mod2_diag <- gamm4(Positive ~ t1_tfinal + s(Age, k=4, bs='cr') +
  s(Age, by=oT1_Tfinal, k=4, bs='cr'), data=clin_df, random=~(1|bblid), REML=TRUE)

all_models_diag <- tab_model(mod1_diag$gam, mod2_diag$gam)

############## Negative ##############
mod1_diag <- gamm4(Negative ~ t1_tfinal + s(Age, k=4, bs='cr'), data=clin_df,
  random=~(1|bblid), REML=TRUE)
mod2_diag <- gamm4(Negative ~ t1_tfinal + s(Age, k=4, bs='cr') +
  s(Age, by=oT1_Tfinal, k=4, bs='cr'), data=clin_df, random=~(1|bblid), REML=TRUE)

all_models_diag <- tab_model(mod1_diag$gam, mod2_diag$gam)

############## Disorganized ##############
mod1_diag <- gamm4(Disorganized ~ t1_tfinal + s(Age, k=4, bs='cr'), data=clin_df,
  random=~(1|bblid), REML=TRUE)
mod2_diag <- gamm4(Disorganized ~ t1_tfinal + s(Age, k=4, bs='cr') +
  s(Age, by=oT1_Tfinal, k=4, bs='cr'), data=clin_df, random=~(1|bblid), REML=TRUE)

all_models_diag <- tab_model(mod1_diag$gam, mod2_diag$gam)

############## General ##############
mod1_diag <- gamm4(General ~ t1_tfinal + s(Age, k=4, bs='cr'), data=clin_df,
  random=~(1|bblid), REML=TRUE)
mod2_diag <- gamm4(General ~ t1_tfinal + s(Age, k=4, bs='cr') +
  s(Age, by=oT1_Tfinal, k=4, bs='cr'), data=clin_df, random=~(1|bblid), REML=TRUE)

all_models_diag <- tab_model(mod1_diag$gam, mod2_diag$gam)

############## GAF_Current ##############
mod1_diag <- gamm4(GAF_Current ~ t1_tfinal + s(Age, k=4, bs='cr'), data=clin_df,
  random=~(1|bblid), REML=TRUE)
mod2_diag <- gamm4(GAF_Current ~ t1_tfinal + s(Age, k=4, bs='cr') +
  s(Age, by=oT1_Tfinal, k=4, bs='cr'), data=clin_df, random=~(1|bblid), REML=TRUE)

all_models_diag <- tab_model(mod1_diag$gam, mod2_diag$gam)

############## GAF_Highest ##############
#mod1_diag <- gamm4(GAF_Highest ~ t1_tfinal + s(Age, k=4, bs='cr'), data=clin_df,
#  random=~(1|bblid), REML=TRUE)
#mod2_diag <- gamm4(GAF_Highest ~ t1_tfinal + s(Age, k=4, bs='cr') +
#  s(Age, by=oT1_Tfinal, k=4, bs='cr'), data=clin_df, random=~(1|bblid), REML=TRUE)

#all_models_diag <- tab_model(mod1_diag$gam, mod2_diag$gam)









#
