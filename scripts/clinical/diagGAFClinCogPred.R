### This script runs models of baseline factor scores and other clinical dimensions
### predicting final PS status
###
### Ellyn Butler
### December 1, 2020 - December 2, 2020

library('dplyr') # Version 1.0.2
library('sjPlot') # Version 2.8.6
library('reshape2') # Version 1.4.4
library('fastDummies') # Version 1.6.3

clin_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/clinical/pnc_longitudinal_diagnosis_n752_202007.csv', stringsAsFactors = TRUE)
demo_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/demographics/baseline/n1601_demographics_go1_20161212.csv', stringsAsFactors = TRUE)
#factor_df <- read.csv('~/Documents/predLongLabels/data/n9498_goassess_itemwise_corrtraits_scores_20161219.csv', stringsAsFactors = TRUE)
ratdiag_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/clinical/n9498_goassess_psych_summary_vars_20131014.csv', stringsAsFactors = TRUE)
cog_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/cognitive/cnb_quickFA_impute_2020-10-20.csv', stringsAsFactors = TRUE)


cog_df <- cog_df[cog_df$Timepoint == 1, ]
names(cog_df)[names(cog_df) == 'EFF_Soln4_MR1'] <- 'SocCog_EFF'
names(cog_df)[names(cog_df) == 'EFF_Soln4_MR2'] <- 'Exec_EFF'
names(cog_df)[names(cog_df) == 'EFF_Soln4_MR3'] <- 'Mem_EFF'
names(cog_df)[names(cog_df) == 'EFF_Soln4_MR4'] <- 'CompCog_EFF'
cog_df <- cog_df[, c('bblid', 'SocCog_EFF', 'Exec_EFF', 'Mem_EFF', 'CompCog_EFF')]

colsforplot <- grep('smry', names(ratdiag_df), value=TRUE)
colsforplot <- colsforplot[!(colsforplot %in% c('smry_prime_pos1', 'smry_prime_tot',
  'smry_prime_pos2', 'smry_psych_overall_rtg', grep('hal', colsforplot, value=TRUE),
  grep('del', colsforplot, value=TRUE)))]

ratdiag_df$ratings_sum <- scale(rowSums(ratdiag_df[, colsforplot]))


final_df <- merge(clin_df, demo_df, by='bblid')
final_df <- merge(final_df, ratdiag_df, by='bblid')
final_df <- merge(final_df, cog_df, by='bblid')
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
final_df$SocCog_EFF <- scale(final_df$SocCog_EFF)
final_df$Exec_EFF <- scale(final_df$Exec_EFF)
final_df$Mem_EFF <- scale(final_df$Mem_EFF)
final_df$CompCog_EFF <- scale(final_df$CompCog_EFF)
final_df$ageAtClinicalAssess1 <- scale(final_df$ageAtClinicalAssess1)
final_df$ageAtClinicalAssess1_sq <- scale((final_df$ageAtClinicalAssess1 - mean(final_df$ageAtClinicalAssess1))^2)


#final_df <- merge(final_df, factor_df)
#final_df[, grep('corr', names(factor_df), value=TRUE)] <- sapply(final_df[, grep('corr', names(factor_df), value=TRUE)], scale)

###############################################################################

################ Exec_EFF
# Sum of ratings, controlling for age
mod0_rat <- glm(PS_final ~ ageAtClinicalAssess1 + ageAtClinicalAssess1_sq + first_diagnosis, family='binomial', data=final_df)
mod1_rat <- glm(PS_final ~ ageAtClinicalAssess1 + ageAtClinicalAssess1_sq + ratings_sum, family='binomial', data=final_df)
mod2_rat <- glm(PS_final ~ ageAtClinicalAssess1 + ageAtClinicalAssess1_sq + Exec_EFF, family='binomial', data=final_df)
mod3_rat <- glm(PS_final ~ ageAtClinicalAssess1 + ageAtClinicalAssess1_sq + first_diagnosis + ratings_sum + Exec_EFF,
  family='binomial', data=final_df)
mod4_rat <- glm(PS_final ~ ageAtClinicalAssess1 + ageAtClinicalAssess1_sq + first_diagnosis + ratings_sum + Exec_EFF + first_diagnosis:ratings_sum + first_diagnosis:Exec_EFF + ratings_sum:Exec_EFF,
  family='binomial', data=final_df)
mod5_rat <- glm(PS_final ~ ageAtClinicalAssess1 + ageAtClinicalAssess1_sq + first_diagnosis*ratings_sum*Exec_EFF,
  family='binomial', data=final_df)

tab_model(mod0_rat, mod1_rat, mod2_rat, mod3_rat, mod4_rat, mod5_rat)


################ CompCog_EFF
# Sum of ratings, controlling for age
mod2_rat <- glm(PS_final ~ ageAtClinicalAssess1 + ageAtClinicalAssess1_sq + CompCog_EFF, family='binomial', data=final_df)
mod3_rat <- glm(PS_final ~ ageAtClinicalAssess1 + ageAtClinicalAssess1_sq + first_diagnosis + ratings_sum + CompCog_EFF,
  family='binomial', data=final_df)
mod4_rat <- glm(PS_final ~ ageAtClinicalAssess1 + ageAtClinicalAssess1_sq + first_diagnosis + ratings_sum + CompCog_EFF + first_diagnosis:ratings_sum + first_diagnosis:CompCog_EFF + ratings_sum:CompCog_EFF,
  family='binomial', data=final_df)
mod5_rat <- glm(PS_final ~ ageAtClinicalAssess1 + ageAtClinicalAssess1_sq + first_diagnosis*ratings_sum*CompCog_EFF,
  family='binomial', data=final_df)

tab_model(mod0_rat, mod1_rat, mod2_rat, mod3_rat, mod4_rat, mod5_rat)

################ Mem_EFF
# Sum of ratings, controlling for age
mod2_rat <- glm(PS_final ~ ageAtClinicalAssess1 + ageAtClinicalAssess1_sq + Mem_EFF, family='binomial', data=final_df)
mod3_rat <- glm(PS_final ~ ageAtClinicalAssess1 + ageAtClinicalAssess1_sq + first_diagnosis + ratings_sum + Mem_EFF,
  family='binomial', data=final_df)
mod4_rat <- glm(PS_final ~ ageAtClinicalAssess1 + ageAtClinicalAssess1_sq + first_diagnosis + ratings_sum + Mem_EFF + first_diagnosis:ratings_sum + first_diagnosis:Mem_EFF + ratings_sum:Mem_EFF,
  family='binomial', data=final_df)
mod5_rat <- glm(PS_final ~ ageAtClinicalAssess1 + ageAtClinicalAssess1_sq + first_diagnosis*ratings_sum*Mem_EFF,
  family='binomial', data=final_df)

tab_model(mod0_rat, mod1_rat, mod2_rat, mod3_rat, mod4_rat, mod5_rat)

################ SocCog_EFF
# Sum of ratings, controlling for age
mod2_rat <- glm(PS_final ~ ageAtClinicalAssess1 + ageAtClinicalAssess1_sq + SocCog_EFF, family='binomial', data=final_df)
mod3_rat <- glm(PS_final ~ ageAtClinicalAssess1 + ageAtClinicalAssess1_sq + first_diagnosis + ratings_sum + SocCog_EFF,
  family='binomial', data=final_df)
mod4_rat <- glm(PS_final ~ ageAtClinicalAssess1 + ageAtClinicalAssess1_sq + first_diagnosis + ratings_sum + SocCog_EFF + first_diagnosis:ratings_sum + first_diagnosis:SocCog_EFF + ratings_sum:SocCog_EFF,
  family='binomial', data=final_df)
mod5_rat <- glm(PS_final ~ ageAtClinicalAssess1 + ageAtClinicalAssess1_sq + first_diagnosis*ratings_sum*SocCog_EFF,
  family='binomial', data=final_df)

tab_model(mod0_rat, mod1_rat, mod2_rat, mod3_rat, mod4_rat, mod5_rat)



#################################### No age ####################################
mod0_rat <- glm(PS_final ~ first_diagnosis, family='binomial', data=final_df)
mod1_rat <- glm(PS_final ~ ratings_sum, family='binomial', data=final_df)
mod2_rat <- glm(PS_final ~ CompCog_EFF, family='binomial', data=final_df)
mod3_rat <- glm(PS_final ~ first_diagnosis + ratings_sum + CompCog_EFF,
  family='binomial', data=final_df)
mod4_rat <- glm(PS_final ~ first_diagnosis + ratings_sum + CompCog_EFF + first_diagnosis:ratings_sum + first_diagnosis:CompCog_EFF + ratings_sum:CompCog_EFF,
  family='binomial', data=final_df)
mod5_rat <- glm(PS_final ~ first_diagnosis*ratings_sum*CompCog_EFF,
  family='binomial', data=final_df)

tab_model(mod0_rat, mod1_rat, mod2_rat, mod3_rat, mod4_rat, mod5_rat)

########################### Age-regressed cognition ###########################


final_df$Exec_EFF_ar <- scale(lm(Exec_EFF ~ ageAtClinicalAssess1 + ageAtClinicalAssess1_sq, data=final_df)$residuals)
mod_exec_ar <- glm(PS_final ~ first_diagnosis + ratings_sum + Exec_EFF_ar,
  family='binomial', data=final_df)

final_df$Mem_EFF_ar <- scale(lm(Mem_EFF ~ ageAtClinicalAssess1 + ageAtClinicalAssess1_sq, data=final_df)$residuals)
mod_mem_ar <- glm(PS_final ~ first_diagnosis + ratings_sum + Mem_EFF_ar,
  family='binomial', data=final_df)

final_df$CompCog_EFF_ar <- scale(lm(CompCog_EFF ~ ageAtClinicalAssess1 + ageAtClinicalAssess1_sq, data=final_df)$residuals)
mod_compcog_ar <- glm(PS_final ~ first_diagnosis + ratings_sum + CompCog_EFF_ar,
  family='binomial', data=final_df)

final_df$SocCog_EFF_ar <- scale(lm(SocCog_EFF ~ ageAtClinicalAssess1 + ageAtClinicalAssess1_sq, data=final_df)$residuals)
mod_soccog_ar <- glm(PS_final ~ first_diagnosis + ratings_sum + SocCog_EFF_ar,
  family='binomial', data=final_df)

tab_model(mod_exec_ar, mod_mem_ar, mod_compcog_ar, mod_soccog_ar)
















  #
