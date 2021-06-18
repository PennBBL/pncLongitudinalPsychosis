### This script runs logistic regressions to predict final PS status using baseline
### number of types of trauma and neighborhood environment.
### Creates Figure 4 and Table 3 in the longitudinal clinical paper.
###
### Ellyn Butler
### September 8, 2020 - September 11, 2020 (changed to envSES on November 30, 2020)
### New diagnostic labels on March 4, 2021

library('dplyr') # Version 1.0.2
library('sjPlot') # Version 2.8.4
library('reshape2') # Version 1.4.4
library('ggplot2') # Version 3.3.2
library('ggpubr') # Version 0.4.0
library('fastDummies') # Version 1.6.1
library('cowplot') # Version 1.1.1

clin_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/clinical/pnc_longitudinal_diagnosis_n749_20210112.csv', stringsAsFactors = TRUE)
demo_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/demographics/baseline/n9498_demographics_go1_20161212.csv', stringsAsFactors = TRUE)
env_df <- read.csv('~/Documents/traumaInformant/data/n9498_go1_environment_factor_scores_tymoore_20150909.csv', stringsAsFactors = TRUE)
trauma_df <- read.csv('~/Documents/traumaInformant/data/PNC_GO1_GOASSESSDataArchiveNontext_DATA_2015-07-14_1157.csv', stringsAsFactors = TRUE)

# Recalculate ntimepoints (a lot of erroneous zeros)
for (bblid in unique(clin_df$bblid)) {
  clin_df[clin_df$bblid == bblid, 'ntimepoints'] <- length(clin_df[clin_df$bblid == bblid, 'ntimepoints'])
}

# Create first/last diagnoses df
clin_df <- clin_df[clin_df$timepoint == 't1' | clin_df$timepoint == paste0('t', clin_df$ntimepoints), ]
clin_df$timepoint <- recode(clin_df$timepoint, 't1'='t1', 't2'='tfinal2',
  't3'='tfinal2', 't4'='tfinal2', 't5'='tfinal2', 't6'='tfinal2')

clin_df$diagnosis <- recode(clin_df$diagnosis, 'psy'='PS')

clin_df <- reshape2::dcast(clin_df, bblid ~ timepoint, value.var='diagnosis')
clin_df$t1_tfinal <- paste(clin_df$t1, clin_df$tfinal2, sep='_')

clin_df$Diagnoses <- recode(clin_df$t1_tfinal, 'TD_TD'='TD-TD', 'TD_other'='TD-OP',
  'TD_PS'='TD-PS', 'other_TD'='OP-TD', 'other_other'='OP-OP', 'other_PS'='OP-PS',
  'PS_TD'='PS-TD', 'PS_other'='PS-OP', 'PS_PS'='PS-PS')
clin_df$Diagnoses <- factor(clin_df$Diagnoses)

trauma_df <- trauma_df[trauma_df$interview_type %in% c('AP', 'MP', 'YPI'),]
trauma_df <- trauma_df[,c('proband_bblid', grep('ptd', names(trauma_df), value=TRUE))]
names(trauma_df)[names(trauma_df) == 'proband_bblid'] <- 'bblid'

final_df <- merge(clin_df, demo_df, by='bblid')
final_df <- merge(final_df, trauma_df, by='bblid')
names(final_df)[names(final_df) == 't1'] <- 'first_diagnosis'
names(final_df)[names(final_df) == 'tfinal2'] <- 'last_diagnosis'
final_df$Female <- recode(final_df$sex, `2`=1, `1`=0)
final_df$Sex <- recode(final_df$sex, `2`='Female', `1`='Male')
final_df$White <- recode(final_df$race, `1`='Yes', .default='No')
final_df$first_diagnosis <- recode(final_df$first_diagnosis, 'other'='OP')
final_df$first_diagnosis <- factor(final_df$first_diagnosis)
final_df$last_diagnosis <- recode(final_df$last_diagnosis, 'other'='OP')
#final_df$t1_tfinal <- recode(final_df$t1_tfinal, 'other_other'='OP_OP',
#  'other_TD'='OP_TD', 'other_PS'='OP_PS', 'TD_other'='TD_OP', 'PS_other'='PS_OP')
final_df <- within(final_df, Diagnoses <- relevel(Diagnoses, ref='TD-TD'))
final_df <- within(final_df, first_diagnosis <- relevel(first_diagnosis, ref='TD'))

ptdvars <- c(paste0('ptd00', 1:4), paste0('ptd00', 6:9))
final_df[, ptdvars] <- sapply(final_df[, ptdvars], na_if, y=9)

final_df <- merge(final_df, env_df)
final_df$PS_final <- recode(final_df$last_diagnosis, 'OP'=0, 'PS'=1, 'TD'=0)

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

#final_df$PS_final_I <- recode(final_df$PS_final, 'Yes'=1, 'No'=0)
trauma_plot <- ggplot(final_df, aes(x=num_type_trauma, y=PS_final,
  group=first_diagnosis, colour=first_diagnosis)) +
  theme_linedraw() + geom_point(shape=1, position=position_jitter(width=.1,
    height=.05)) + ylab('PS at Final Assessment') + xlab('Number of Types of Trauma') +
  labs(colour='First Diagnosis') +
  scale_colour_manual(values = c('green3', 'goldenrod2', 'red')) +
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
env_plot <- ggplot(final_df, aes(x=envSES, y=PS_final,
  group=first_diagnosis, colour=first_diagnosis)) +
  theme_linedraw() + geom_point(shape=1, position=position_jitter(width=.1,
    height=.05)) + ylab('PS at Final Assessment') + xlab('Neighborhood Environment') +
  labs(colour='First Diagnosis') +
  scale_colour_manual(values = c('green3', 'goldenrod2', 'red')) +
  stat_smooth(method='glm', method.args=list(family='binomial'), se=FALSE) +
  theme(legend.position='bottom')

pdf(file='~/Documents/pncLongitudinalPsychosis/plots/envLogistic.pdf', width=6, height=5)
env_plot
dev.off()

# Built figure
diag_legend <- get_legend(env_plot)

trauma_plot <- trauma_plot + theme(legend.position='none')
env_plot <- env_plot + theme(legend.position='none')
grid_plot <- cowplot::plot_grid(
  cowplot::plot_grid(trauma_plot, env_plot, labels=c('A', 'B')),
  diag_legend, rel_heights=c(4, 1), nrow=2, ncol=1)

pdf(file='~/Documents/pncLongitudinalPsychosis/plots/env_grid_paper.pdf', width=7, height=4)
print(grid_plot)
dev.off()


###############################################################################
###############################################################################
###############################################################################

###### Does the interaction between trauma and environment predict? ######

mod2.3 <- glm(PS_final ~ first_diagnosis + num_type_trauma + envSES, family='binomial', data=final_df)

mod3.3 <- glm(PS_final ~ first_diagnosis + num_type_trauma*envSES, family='binomial', data=final_df)

mod4.3 <- glm(PS_final ~ first_diagnosis*num_type_trauma*envSES, family='binomial', data=final_df)

mod5.3 <- glm(PS_final ~ first_diagnosis*num_type_trauma*envSES*sex, family='binomial', data=final_df) #Wildly overfit

print(tab_model(mod1, mod2.2, mod2, mod2.3, mod3.3, mod4.3, p.adjust='fdr',
  file='~/Documents/pncLongitudinalPsychosis/results/table_prediction.html'))

###############################################################################

final_df$Age <- final_df$ageAtClinicalAssess1

# Sensitivity Table
mod1b <- glm(PS_final ~ Age + Sex + White + first_diagnosis, family='binomial', data=final_df)

mod2.2b <- glm(PS_final ~ Age + Sex + White + first_diagnosis + envSES, family='binomial', data=final_df)

mod2b <- glm(PS_final ~ Age + Sex + White + first_diagnosis + num_type_trauma, family='binomial', data=final_df)

mod2.3b <- glm(PS_final ~ Age + Sex + White + first_diagnosis + num_type_trauma + envSES, family='binomial', data=final_df)

mod3.3b <- glm(PS_final ~ Age + Sex + White + first_diagnosis + num_type_trauma*envSES, family='binomial', data=final_df)

mod4.3b <- glm(PS_final ~ Age + Sex + White + first_diagnosis*num_type_trauma*envSES, family='binomial', data=final_df)

print(tab_model(mod1b, mod2.2b, mod2b, mod2.3b, mod3.3b, mod4.3b, p.adjust='fdr',
  file='~/Documents/pncLongitudinalPsychosis/results/table_prediction_sensitivity.html'))












    #
