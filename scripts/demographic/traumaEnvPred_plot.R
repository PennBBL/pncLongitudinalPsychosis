### This script makes plot of baseline demographic and trauma features by
### 3x3 trajectories
###
### Ellyn Butler
### September 8, 2020 - September 11, 2020 (changed to envSES on November 30, 2020)

library('dplyr') # Version 1.0.2
library('sjPlot') # Version 2.8.4
library('reshape2') # Version 1.4.4
library('ggplot2') # Version 3.3.2
library('ggpubr') # Version 0.4.0
library('fastDummies') # Version 1.6.1

clin_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/clinical/pnc_longitudinal_diagnosis_n752_202007.csv')
demo_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/demographics/baseline/n1601_demographics_go1_20161212.csv')
env_df <- read.csv('~/Documents/traumaInformant/data/n9498_go1_environment_factor_scores_tymoore_20150909.csv')
trauma_df <- read.csv('~/Documents/traumaInformant/data/PNC_GO1_GOASSESSDataArchiveNontext_DATA_2015-07-14_1157.csv')

trauma_df <- trauma_df[trauma_df$interview_type %in% c('AP', 'MP', 'YPI'),]
trauma_df <- trauma_df[,c('proband_bblid', grep('ptd', names(trauma_df), value=TRUE))]
names(trauma_df)[names(trauma_df) == 'proband_bblid'] <- 'bblid'

final_df <- merge(clin_df, demo_df, by='bblid')
final_df <- merge(final_df, trauma_df, by='bblid')
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

ptdvars <- c(paste0('ptd00', 1:4), paste0('ptd00', 6:9))
final_df[, ptdvars] <- sapply(final_df[, ptdvars], na_if, y=9)

final_df <- merge(final_df, env_df)
final_df$PS_final <- recode(final_df$tfinal, 'other'='No', 'PS'='Yes',
  'psy'='Yes', 'TD'='No')

final_df$num_type_trauma <- rowSums(final_df[, ptdvars])

######### Does # of trauma types predict final time point PS status? #########

mod1 <- glm(PS_final ~ first_diagnosis, family='binomial', data=final_df)

mod2 <- glm(PS_final ~ first_diagnosis + num_type_trauma, family='binomial', data=final_df)

#Answer: Yes

# Does baseline diagnosis moderate the relationship between final time
# point PS status and number of types of traumas such that there is a stronger
# relationship between number of trauma types and final PS status among those
# who were PS at baseline than those who were OP or TD?

mod3 <- glm(PS_final ~ num_type_trauma*first_diagnosis, family='binomial', data=final_df)

#Answer: No

all_models_trauma <- tab_model(mod1, mod2, mod3)

final_df$PS_final_I <- recode(final_df$PS_final, 'Yes'=1, 'No'=0)
trauma_plot <- ggplot(final_df, aes(x=num_type_trauma, y=PS_final_I,
  group=first_diagnosis, colour=first_diagnosis)) +
  theme_linedraw() + geom_point(shape=1, position=position_jitter(width=.1,
    height=.05)) + ylab('PS final') +
  scale_colour_manual(values = c('slategray2', 'deeppink3', 'black')) +
  stat_smooth(method='glm', method.args=list(family='binomial'), se=FALSE) +
  theme(legend.position='bottom')

pdf(file='~/Documents/pncLongitudinalPsychosis/plots/traumaLogistic.pdf', width=6, height=5)
trauma_plot
dev.off()

###### Does neighborhood environment predict final time point PS status? ######


mod2.2 <- glm(PS_final ~ first_diagnosis + envSES, family='binomial', data=final_df)

#Answer: No

# Does baseline diagnosis moderate the relationship between final time
# point PS status and number of types of traumas such that there is a stronger
# relationship between number of trauma types and final PS status among those
# who were PS at baseline than those who were OP or TD?

mod3.2 <- glm(PS_final ~ envSES*first_diagnosis, family='binomial', data=final_df)

#Answer: No

all_models_env <- tab_model(mod1, mod2.2, mod3.2)

env_plot <- ggplot(final_df, aes(x=envSES, y=PS_final_I,
  group=first_diagnosis, colour=first_diagnosis)) +
  theme_linedraw() + geom_point(shape=1, position=position_jitter(width=.1,
    height=.05)) + ylab('PS final') +
  scale_colour_manual(values = c('slategray2', 'deeppink3', 'black')) +
  stat_smooth(method='glm', method.args=list(family='binomial'), se=FALSE) +
  theme(legend.position='bottom')

pdf(file='~/Documents/pncLongitudinalPsychosis/plots/envLogistic.pdf', width=6, height=5)
env_plot
dev.off()




###############################################################################
###############################################################################
###############################################################################

###### Does the interaction between trauma and environment predict? ######

mod2.3 <- glm(PS_final ~ first_diagnosis + num_type_trauma + envSES, family='binomial', data=final_df)

mod3.3 <- glm(PS_final ~ first_diagnosis + num_type_trauma*envSES, family='binomial', data=final_df)

mod4.3 <- glm(PS_final ~ first_diagnosis*num_type_trauma*envSES, family='binomial', data=final_df)

mod5.3 <- glm(PS_final ~ first_diagnosis*num_type_trauma*envSES*sex, family='binomial', data=final_df) #Wildly overfit

all_models_both <- tab_model(mod1, mod2, mod2.3, mod3.3, mod4.3)

##### Within first diagnosis groups

#### TD
for (diag in c('TD', 'OP', 'PS')) {
  mod1 <- glm(PS_final ~ sex, family='binomial',
    data=final_df[final_df$first_diagnosis == diag, ])
  mod2 <- glm(PS_final ~ sex + num_type_trauma,
    family='binomial', data=final_df[final_df$first_diagnosis == diag, ])
  mod3 <- glm(PS_final ~ sex + num_type_trauma + envSES,
    family='binomial', data=final_df[final_df$first_diagnosis == diag, ])
  mod4 <- glm(PS_final ~ sex + num_type_trauma + envSES +
    sex:num_type_trauma + num_type_trauma:envSES + sex:envSES,
    family='binomial', data=final_df[final_df$first_diagnosis == diag, ])
  mod5 <- glm(PS_final ~ sex*num_type_trauma*envSES,
      family='binomial', data=final_df[final_df$first_diagnosis == diag, ])

  mods <- tab_model(mod1, mod2, mod3, mod4, mod5)
  assign(paste0('mods_', diag), mods)
}
