### This script runs models of baseline factor scores and other clinical dimensions
### predicting final PS status
###
### Ellyn Butler
### November 5, 2020 - November 30, 2020

library('dplyr') # Version 1.0.2
library('sjPlot') # Version 2.8.4
library('reshape2') # Version 1.4.4
library('fastDummies') # Version 1.6.1

clin_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/clinical/pnc_longitudinal_diagnosis_n752_202007.csv')
demo_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/demographics/baseline/n1601_demographics_go1_20161212.csv')
factor_df <- read.csv('~/Documents/predLongLabels/data/n9498_goassess_itemwise_corrtraits_scores_20161219.csv')
ratdiag_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/clinical/n9498_goassess_psych_summary_vars_20131014.csv')

colsforplot <- grep('smry', names(ratdiag_df), value=TRUE)
colsforplot <- colsforplot[!(colsforplot %in% c('smry_prime_pos1', 'smry_prime_tot',
  'smry_prime_pos2', 'smry_psych_overall_rtg', grep('hal', colsforplot, value=TRUE),
  grep('del', colsforplot, value=TRUE)))]

ratdiag_df$ratings_sum <- scale(rowSums(ratdiag_df[, colsforplot]))


final_df <- merge(clin_df, demo_df, by='bblid')
final_df <- merge(final_df, ratdiag_df, by='bblid')
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

final_df$PS_final <- recode(final_df$tfinal, 'other'='No', 'PS'='Yes',
  'psy'='Yes', 'TD'='No')

final_df$sex <- recode(final_df$sex, `2`='Female', `1`='Male')
final_df$sex <- as.factor(final_df$sex)
final_df <- within(final_df, sex <- relevel(sex, ref='Male'))

final_df <- merge(final_df, factor_df)
final_df[, grep('corr', names(factor_df), value=TRUE)] <- sapply(final_df[, grep('corr', names(factor_df), value=TRUE)], scale)

###############################################################################

# Correlated traits
mod0 <- glm(PS_final ~ first_diagnosis + mood_corrtraitsv2, family='binomial', data=final_df)
mod1 <- glm(PS_final ~ first_diagnosis + psychosis_corrtraitsv2, family='binomial', data=final_df)
mod2 <- glm(PS_final ~ first_diagnosis + externalizing_corrtraitsv2, family='binomial', data=final_df)
mod3 <- glm(PS_final ~ first_diagnosis + fear_corrtraitsv2, family='binomial', data=final_df)
mod4 <- glm(PS_final ~ first_diagnosis + mood_corrtraitsv2 + psychosis_corrtraitsv2 + externalizing_corrtraitsv2 + fear_corrtraitsv2, family='binomial', data=final_df)


# Sum of ratings
mod0_rat <- glm(PS_final ~ first_diagnosis, family='binomial', data=final_df)
mod1_rat <- glm(PS_final ~ ratings_sum, family='binomial', data=final_df)
mod2_rat <- glm(PS_final ~ first_diagnosis + ratings_sum, family='binomial', data=final_df)
tab_model(mod0_rat, mod1_rat, mod2_rat)

tab_model(mod0_rat, mod0, mod1, mod2, mod3, mod4)

# Indicators for psychosis, externalizing, internalizing, other
