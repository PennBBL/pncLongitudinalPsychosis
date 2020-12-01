### This script runs GAMMs to see if there are main effects and interactions
### for sex and longitudinal diagnostic labels for the social cognition factor
### score.
###
### Ellyn Butler
### September 8, 2020 - October 21, 2020 (new data October 20, 2020)

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

cnb_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/cognitive/cnb_quickFA_impute_2020-10-20.csv')

names(cnb_df)[names(cnb_df) == 'EFF_Soln4_MR1'] <- 'SocCog_EFF'
names(cnb_df)[names(cnb_df) == 'EFF_Soln4_MR2'] <- 'Exec_EFF'
names(cnb_df)[names(cnb_df) == 'EFF_Soln4_MR3'] <- 'Mem_EFF'
names(cnb_df)[names(cnb_df) == 'EFF_Soln4_MR4'] <- 'CompCog_EFF'

cnb_df$sex <- relevel(cnb_df$sex, 'Male')
cnb_df$t1_tfinal <- relevel(cnb_df$t1_tfinal, 'TD_TD')

cnb_df$oSex <- ordered(cnb_df$sex, c('Male', 'Female'))

cnb_df$oT1_Tfinal <- ordered(cnb_df$t1_tfinal, c('TD_TD', 'OP_OP', 'OP_PS',
  'OP_TD', 'PS_OP', 'PS_PS', 'PS_TD', 'TD_OP', 'TD_PS'))

cnb_df$sex_t1_tfinal <- paste(cnb_df$sex, cnb_df$t1_tfinal, sep='_')
cnb_df$oSex_oT1_Tfinal <- ordered(cnb_df$sex_t1_tfinal, c('Male_TD_TD',
  'Male_OP_OP', 'Male_OP_PS', 'Male_OP_TD', 'Male_PS_OP', 'Male_PS_PS',
  'Male_PS_TD', 'Male_TD_OP', 'Male_TD_PS', 'Female_TD_TD', 'Female_OP_OP',
  'Female_OP_PS', 'Female_OP_TD', 'Female_PS_OP', 'Female_PS_PS', 'Female_PS_TD',
  'Female_TD_OP', 'Female_TD_PS'))


################################## GAMMs ##################################

############## Social Cognition Efficiency ##############
#### Sex Models
mod1_sex <- gamm4(SocCog_EFF ~ sex + s(Age, k=4, bs='cr'), data=cnb_df,
  random=~(1|bblid), REML=TRUE)
mod2_sex <- gamm4(SocCog_EFF ~ sex + s(Age, k=4, bs='cr') +
  s(Age, by=oSex, k=4, bs='cr'), data=cnb_df, random=~(1|bblid), REML=TRUE)

all_models_sex <- tab_model(mod1_sex$gam, mod2_sex$gam)

#### Diagnosis Models
mod1_diag <- gamm4(SocCog_EFF ~ t1_tfinal + s(Age, k=4, bs='cr'), data=cnb_df,
  random=~(1|bblid), REML=TRUE)
mod2_diag <- gamm4(SocCog_EFF ~ t1_tfinal + s(Age, k=4, bs='cr') +
  s(Age, by=oT1_Tfinal, k=4, bs='cr'), data=cnb_df, random=~(1|bblid), REML=TRUE)

all_models_diag <- tab_model(mod1_diag$gam, mod2_diag$gam)

#### Sex by Diagnosis Models
mod1_sex_diag <- gamm4(SocCog_EFF ~ sex + t1_tfinal + s(Age, k=4, bs='cr') +
  s(Age, by=oSex, k=4, bs='cr'), data=cnb_df, random=~(1|bblid), REML=TRUE)

mod2_sex_diag <- gamm4(SocCog_EFF ~ sex + t1_tfinal + s(Age, k=4, bs='cr') +
  s(Age, by=oSex, k=4, bs='cr') + s(Age, by=oT1_Tfinal, k=4, bs='cr'),
  data=cnb_df, random=~(1|bblid), REML=TRUE)

mod3_sex_diag <- gamm4(SocCog_EFF ~ sex*t1_tfinal + s(Age, k=4, bs='cr') +
  s(Age, by=oSex_oT1_Tfinal, k=4, bs='cr'), data=cnb_df, random=~(1|bblid), REML=TRUE)

all_models_sex_diag <- tab_model(mod1_sex$gam, mod2_sex$gam, mod1_sex_diag$gam,
  mod2_sex_diag$gam, mod3_sex_diag$gam)


############## Executive Efficiency ##############
#### Sex Models
mod1_sex <- gamm4(Exec_EFF ~ sex + s(Age, k=4, bs='cr'), data=cnb_df,
  random=~(1|bblid), REML=TRUE)
mod2_sex <- gamm4(Exec_EFF ~ sex + s(Age, k=4, bs='cr') +
  s(Age, by=oSex, k=4, bs='cr'), data=cnb_df, random=~(1|bblid), REML=TRUE)

all_models_sex <- tab_model(mod1_sex$gam, mod2_sex$gam)

#### Diagnosis Models
mod1_diag <- gamm4(Exec_EFF ~ t1_tfinal + s(Age, k=4, bs='cr'), data=cnb_df,
  random=~(1|bblid), REML=TRUE)
mod2_diag <- gamm4(Exec_EFF ~ t1_tfinal + s(Age, k=4, bs='cr') +
  s(Age, by=oT1_Tfinal, k=4, bs='cr'), data=cnb_df, random=~(1|bblid), REML=TRUE)

mod3_sex_diag <- gamm4(Exec_EFF ~ sex*t1_tfinal + s(Age, k=4, bs='cr') +
  s(Age, by=oSex_oT1_Tfinal, k=4, bs='cr'), data=cnb_df, random=~(1|bblid), REML=TRUE)

all_models_diag <- tab_model(mod1_diag$gam, mod2_diag$gam)


############## Memory Efficiency ##############
#### Sex Models
mod1_sex <- gamm4(Mem_EFF ~ sex + s(Age, k=4, bs='cr'), data=cnb_df,
  random=~(1|bblid), REML=TRUE)
mod2_sex <- gamm4(Mem_EFF ~ sex + s(Age, k=4, bs='cr') +
  s(Age, by=oSex, k=4, bs='cr'), data=cnb_df, random=~(1|bblid), REML=TRUE)

all_models_sex <- tab_model(mod1_sex$gam, mod2_sex$gam)

#### Diagnosis Models
mod1_diag <- gamm4(Mem_EFF ~ t1_tfinal + s(Age, k=4, bs='cr'), data=cnb_df,
  random=~(1|bblid), REML=TRUE)
mod2_diag <- gamm4(Mem_EFF ~ t1_tfinal + s(Age, k=4, bs='cr') +
  s(Age, by=oT1_Tfinal, k=4, bs='cr'), data=cnb_df, random=~(1|bblid), REML=TRUE)

mod3_sex_diag <- gamm4(Mem_EFF ~ sex*t1_tfinal + s(Age, k=4, bs='cr') +
  s(Age, by=oSex_oT1_Tfinal, k=4, bs='cr'), data=cnb_df, random=~(1|bblid), REML=TRUE)

all_models_diag <- tab_model(mod1_diag$gam, mod2_diag$gam)


############## Complex Cognition Efficiency ##############
#### Sex Models
mod1_sex <- gamm4(CompCog_EFF ~ sex + s(Age, k=4, bs='cr'), data=cnb_df,
  random=~(1|bblid), REML=TRUE)
mod2_sex <- gamm4(CompCog_EFF ~ sex + s(Age, k=4, bs='cr') +
  s(Age, by=oSex, k=4, bs='cr'), data=cnb_df, random=~(1|bblid), REML=TRUE)

all_models_sex <- tab_model(mod1_sex$gam, mod2_sex$gam)

#### Diagnosis Models
mod1_diag <- gamm4(CompCog_EFF ~ t1_tfinal + s(Age, k=4, bs='cr'), data=cnb_df,
  random=~(1|bblid), REML=TRUE)
mod2_diag <- gamm4(CompCog_EFF ~ t1_tfinal + s(Age, k=4, bs='cr') +
  s(Age, by=oT1_Tfinal, k=4, bs='cr'), data=cnb_df, random=~(1|bblid), REML=TRUE)

mod3_sex_diag <- gamm4(CompCog_EFF ~ sex*t1_tfinal + s(Age, k=4, bs='cr') +
  s(Age, by=oSex_oT1_Tfinal, k=4, bs='cr'), data=cnb_df, random=~(1|bblid), REML=TRUE)

all_models_diag <- tab_model(mod1_diag$gam, mod2_diag$gam)






#
