### This script plots just the fits for the SIPS dimensions
###
### Ellyn Butler
### August 27, 2020 - March 4, 2021

set.seed(20)

library('dplyr') # V 1.0.2
library('reshape2') # V 1.4.4
library('ggplot2') # V 3.3.2
library('ggpubr') # V 0.4.0
library('lme4') # V 1.1-26
library('gamm4') # V 0.2-6
library('stringr') # V 1.4.0
library('MASS') # V 7.3-53
library('broom') # V 0.7.2
library('sjPlot') # V 2.8.6
library('cowplot') # V 1.1.1
library('pbkrtest') # V 0.5-0.1
#library('gratia') # V 0.4.0

# Declare functions
getAssessmentNumber <- function(i, dataf) {
  bblid <- dataf[i, 'bblid']
  ages <- dataf[dataf$bblid == bblid, 'age']
  ages <- sort(unique(ages))
  which(ages == dataf[i, 'age'])
}

getUpperLowerCI <- function(i) {
  sorted_vec <- unname(sort(fits[i,]))
  lower <- sorted_vec[round(.025*length(sorted_vec))]
  upper <- sorted_vec[round(.975*length(sorted_vec))]
  c(lower, upper)
}

# Read in data
clin_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/clinical/clean_sips_nodup.csv')
diag_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/clinical/pnc_longitudinal_diagnosis_n749_20210112.csv')

# Recalculate ntimepoints (a lot of erroneous zeros)
for (bblid in unique(diag_df$bblid)) {
  diag_df[diag_df$bblid == bblid, 'ntimepoints'] <- length(diag_df[diag_df$bblid == bblid, 'ntimepoints'])
}

# Create first/last diagnoses df
diag_df <- diag_df[diag_df$timepoint == 't1' | diag_df$timepoint == paste0('t', diag_df$ntimepoints), ]
diag_df$timepoint <- recode(diag_df$timepoint, 't1'='t1', 't2'='tfinal2',
  't3'='tfinal2', 't4'='tfinal2', 't5'='tfinal2', 't6'='tfinal2')

diag_df$diagnosis <- recode(diag_df$diagnosis, 'psy'='PS')

diag_df <- reshape2::dcast(diag_df, bblid ~ timepoint, value.var='diagnosis')
diag_df$t1_tfinal <- paste(diag_df$t1, diag_df$tfinal2, sep='_')

diag_df$Diagnoses <- recode(diag_df$t1_tfinal, 'TD_TD'='TD-TD', 'TD_other'='TD-OP',
  'TD_PS'='TD-PS', 'other_TD'='OP-TD', 'other_other'='OP-OP', 'other_PS'='OP-PS',
  'PS_TD'='PS-TD', 'PS_other'='PS-OP', 'PS_PS'='PS-PS')

# Name/relevel variables as necessary
diag_df$Diagnoses <- factor(diag_df$Diagnoses)
diag_df$Diagnoses <- relevel(diag_df$Diagnoses, 'TD-TD')

diag_df$oDiagnoses <- ordered(diag_df$Diagnoses, c('TD-TD', 'OP-OP', 'OP-PS',
  'OP-TD', 'PS-OP', 'PS-PS', 'PS-TD', 'TD-OP', 'TD-PS'))

clin_df <- merge(clin_df, diag_df, by='bblid')

clin_df$White <- recode(clin_df$race, `1`='Yes', `2`='No', `3`='No', `4`='No', `5`='No')

names(clin_df)[names(clin_df) == 'gaf_c'] <- 'Global'
names(clin_df)[names(clin_df) == 'age'] <- 'Age'
names(clin_df)[names(clin_df) == 'sex'] <- 'Sex'

############################ Plot SIPS & functioning ############################

plotcols <- c('Positive', 'Negative', 'Disorganized', 'General', 'Global') #, 'Social', 'Role')

diags <- as.character(unique(clin_df$Diagnoses)[!(unique(clin_df$Diagnoses) %in% c('TD-TD', 'PS-PS'))])

i=1
for (score in plotcols) {
  clin_score_df <- clin_df
  clin_score_df <- clin_score_df[!is.na(clin_score_df[, score]), ]
  row.names(clin_score_df) <- 1:nrow(clin_score_df)

  mod1b <- gamm4(as.formula(paste(score, "~ Diagnoses + s(Age, k=4, bs='cr') +
    s(Age, by=oDiagnoses, k=10, bs='cr')")), data=clin_score_df, random=~(1|bblid), REML=TRUE)

  mod2b <- gamm4(as.formula(paste(score, "~ Sex + White + Diagnoses + s(Age, k=4, bs='cr') +
    s(Age, by=oDiagnoses, k=10, bs='cr')")), data=clin_score_df, random=~(1|bblid), REML=TRUE)

  print(tab_model(mod1b$gam, file=paste0('~/Documents/pncLongitudinalPsychosis/results/table_', score, '.html')))
  assign(paste0(score, '_model'), mod1b$gam)
  print(tab_model(mod2b$gam, file=paste0('~/Documents/pncLongitudinalPsychosis/results/tableSensitivity_', score, '.html')))
  assign(paste0(score, 'Sensitivity_model'), mod2b$gam)

  lp <- predict(mod1b$gam, newdata=clin_score_df, type='lpmatrix')
  coefs <- coef(mod1b$gam)
  vc <- vcov(mod1b$gam)

  sim <- mvrnorm(1000, mu = coefs, Sigma = vc)
  fits <- lp %*% t(sim)

  cis <- t(sapply(1:nrow(fits), getUpperLowerCI))
  cis <- data.frame(cis)
  names(cis) <- c('LCI', 'UCI')
  clin_score_df <- cbind(clin_score_df, cis)

  clin_score_df$predgamm <- predict(mod1b$gam)

  for (group in diags) {
    subtit <- paste0('Sessions: TD-TD=', nrow(clin_score_df[clin_score_df$Diagnoses == 'TD-TD',]),
      ', ', group, '=', nrow(clin_score_df[clin_score_df$Diagnoses == group,]),
      ', PS-PS=', nrow(clin_score_df[clin_score_df$Diagnoses == 'PS-PS',]))
    group_df <- clin_score_df[clin_score_df$Diagnoses %in% c('TD-TD', group, 'PS-PS'), ]
    group_df$Diagnoses <- ordered(group_df$Diagnoses, c('TD-TD', group, 'PS-PS'))

    intervals_df <- reshape2::melt(group_df, c('bblid', 'Age', 'Diagnoses'), c('LCI', 'UCI'))
    clin_plot <- ggplot(group_df, aes_string(x='Age', y=score, color='Diagnoses')) +
      theme_linedraw() + xlim(10, 30) +
      labs(title=score, subtitle=subtit, y='Score (95% CI)', color='Diagnoses') +
      scale_color_manual(values=c('green3', 'goldenrod2', 'red')) +
      theme(legend.position = 'bottom', plot.title=element_text(size=14, face="bold"),
        plot.subtitle=element_text(size=8)) +
      geom_line(aes(y=predgamm), size=1) +
      geom_ribbon(data=clin_score_df[clin_score_df$Diagnoses == 'TD-TD',],
        aes(ymin=LCI, ymax=UCI), alpha=.2, fill='green3', show.legend=FALSE) +
      geom_ribbon(data=clin_score_df[clin_score_df$Diagnoses == group,],
        aes(ymin=LCI, ymax=UCI), alpha=.2, fill='goldenrod2', show.legend=FALSE) +
      geom_ribbon(data=clin_score_df[clin_score_df$Diagnoses == 'PS-PS',],
        aes(ymin=LCI, ymax=UCI), alpha=.2, fill='red', show.legend=FALSE)

    assign(paste0(score, '_', gsub('-', '_', group), '_plot'), clin_plot)

    pdf(file=paste0('~/Documents/pncLongitudinalPsychosis/plots/clin_', score, '_', group, '_clean.pdf'), width=4, height=4)
    print(clin_plot)
    dev.off()
  }
  i=i+1
}

# Add statistics to figures (or maybe not, just create a mega table)
print(tab_model(Positive_model, Negative_model, Disorganized_model, General_model,
  Global_model, p.adjust='fdr', file='~/Documents/pncLongitudinalPsychosis/results/clin_mega.html'))
print(tab_model(PositiveSensitivity_model, NegativeSensitivity_model,
  DisorganizedSensitivity_model, GeneralSensitivity_model, GlobalSensitivity_model,
  p.adjust='fdr', file='~/Documents/pncLongitudinalPsychosis/results/clinSensitivity_mega.html'))

# Build 2x2 with legend (SIPS)
diag_legend <- get_legend(Positive_OP_OP_plot)

Positive_OP_OP_plot <- Positive_OP_OP_plot + theme(legend.position='none')
Negative_OP_OP_plot <- Negative_OP_OP_plot + theme(legend.position='none')
Disorganized_OP_OP_plot <- Disorganized_OP_OP_plot + theme(legend.position='none')
General_OP_OP_plot <- General_OP_OP_plot + theme(legend.position='none')

grid_plot <- cowplot::plot_grid(
  cowplot::plot_grid(Positive_OP_OP_plot, Negative_OP_OP_plot, labels=c('A', 'B')),
  cowplot::plot_grid(Disorganized_OP_OP_plot, General_OP_OP_plot, labels=c('C', 'D')),
  diag_legend, rel_heights=c(4, 4, 1), nrow=3, ncol=1)

pdf(file='~/Documents/pncLongitudinalPsychosis/plots/sips_grid_paper.pdf', width=7, height=8)
print(grid_plot)
dev.off()



# Build 1x3 with legend (functioning) - NOT ENOUGH CORNBLATT DATA
#Global_OP_OP_plot <- Global_OP_OP_plot + theme(legend.position='none')
#Social_OP_OP_plot <- Social_OP_OP_plot + theme(legend.position='none')
#Role_OP_OP_plot <- Role_OP_OP_plot + theme(legend.position='none')

#grid_plot <- cowplot::plot_grid(
#  cowplot::plot_grid(Global_OP_OP_plot, Social_OP_OP_plot, Role_OP_OP_plot,
#  labels=c('A', 'B', 'C'), ncol=3), diag_legend, rel_heights=c(4, 1), nrow=2, ncol=1)

#pdf(file='~/Documents/pncLongitudinalPsychosis/plots/func_grid_paper.pdf', width=10.5, height=4.444)
#print(grid_plot)
#dev.off()






















#
