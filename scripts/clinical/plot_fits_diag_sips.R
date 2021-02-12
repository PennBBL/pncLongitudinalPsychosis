### This script plots just the fits for the SIPS dimensions
###
### Ellyn Butler
### August 27, 2020 - February 10, 2021

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

clin_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/clinical/clean_sips.csv')
diag_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/clinical/pnc_longitudinal_diagnosis_n752_202007.csv')
diag_df$Diagnoses <- recode(diag_df$t1_tfinal, 'TD_TD'='TD-TD', 'TD_other'='TD-OP',
  'TD_PS'='TD-PS', 'other_TD'='OP-TD', 'other_other'='OP-OP', 'other_PS'='OP-PS',
  'PS_TD'='PS-TD', 'PS_other'='PS-OP', 'PS_PS'='PS-PS')

final_df <- merge(clin_df, diag_df, by='bblid')

final_df$Diagnoses <- factor(final_df$Diagnoses)
final_df$Diagnoses <- relevel(final_df$Diagnoses, 'TD-TD')

final_df$oDiagnoses <- ordered(final_df$Diagnoses, c('TD-TD', 'OP-OP', 'OP-PS',
  'OP-TD', 'PS-OP', 'PS-PS', 'PS-TD', 'TD-OP', 'TD-PS'))

names(final_df)[names(final_df) == 'gaf_c'] <- 'GAF'
#names(final_df)[names(final_df) == 'gaf_h'] <- 'GAF_Highest'
names(final_df)[names(final_df) == 'age'] <- 'Age'

getUpperLowerCI <- function(i) {
  sorted_vec <- unname(sort(fits[i,]))
  lower <- sorted_vec[round(.025*length(sorted_vec))]
  upper <- sorted_vec[round(.975*length(sorted_vec))]
  c(lower, upper)
}

############################# Plot SIPS Dimensions #############################

plotcols <- c('Positive', 'Negative', 'Disorganized', 'General', 'GAF')

diags <- as.character(unique(final_df$Diagnoses)[!(unique(final_df$Diagnoses) %in% c('TD-TD', 'PS-PS'))])

for (score in plotcols) {
  clin_score_df <- final_df
  clin_score_df <- clin_score_df[!is.na(clin_score_df[, score]), ]
  row.names(clin_score_df) <- 1:nrow(clin_score_df)

  mod1b <- gamm4(as.formula(paste(score, "~ Diagnoses + s(Age, k=4, bs='cr') +
    s(Age, by=oDiagnoses, k=10, bs='cr')")), data=clin_score_df, random=~(1|bblid), REML=TRUE)

  print(tab_model(mod1b$gam, file=paste0('~/Documents/pncLongitudinalPsychosis/results/table_', score, '.html')))
  assign(paste0(score, '_model'), mod1b$gam)

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

      if (score == 'GAF') {
        clin_plot <- clin_plot + labs(title='Global Assessment of Functioning', subtitle=subtit, y='Score (95% CI)', color='Diagnoses')
      } else {
        clin_plot <- clin_plot + labs(title=score, subtitle=subtit, y='Score (95% CI)', color='Diagnoses')
      }

    assign(paste0(score, '_', gsub('-', '_', group), '_plot'), clin_plot)

    pdf(file=paste0('~/Documents/pncLongitudinalPsychosis/plots/clin_', score, '_', group, '_clean.pdf'), width=4, height=4)
    print(clin_plot)
    dev.off()
  }
}

# Add statistics to figures (or maybe not, just create a mega table)
print(tab_model(Positive_model, Negative_model, Disorganized_model, General_model, GAF_model,
  p.adjust='fdr', file='~/Documents/pncLongitudinalPsychosis/results/clin_mega.html'))


# Build 2x2 with legend
diag_legend <- get_legend(Positive_OP_OP_plot)

Positive_OP_OP_plot <- Positive_OP_OP_plot + theme(legend.position='none')
Negative_OP_OP_plot <- Negative_OP_OP_plot + theme(legend.position='none')
Disorganized_OP_OP_plot <- Disorganized_OP_OP_plot + theme(legend.position='none')
General_OP_OP_plot <- General_OP_OP_plot + theme(legend.position='none')

grid_plot <- cowplot::plot_grid(
  cowplot::plot_grid(Positive_OP_OP_plot, Negative_OP_OP_plot, labels=c('A', 'B')),
  cowplot::plot_grid(Disorganized_OP_OP_plot, General_OP_OP_plot, labels=c('C', 'D')),
  diag_legend, rel_heights=c(4, 4, 1), nrow=3, ncol=1)

pdf(file='~/Documents/pncLongitudinalPsychosis/plots/clin_grid_paper.pdf', width=7, height=8)
print(grid_plot)
dev.off()

























#
