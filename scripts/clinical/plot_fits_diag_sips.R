### This script plots just the fits for the SIPS dimensions
###
### Ellyn Butler
### August 27, 2020 - October 20, 2020

set.seed(20)

library('dplyr')
library('reshape2')
library('ggplot2')
library('ggpubr')
library('lme4')
library('gamm4')
library('stringr')
library('MASS')
library('broom')
library('sjPlot')

clin_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/clinical/clean_sips.csv')
diag_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/clinical/pnc_longitudinal_diagnosis_n752_202007.csv')
diag_df$t1_tfinal <- recode(diag_df$t1_tfinal, 'TD_TD'='TD_TD', 'TD_other'='TD_OP',
  'TD_PS'='TD_PS', 'other_TD'='OP_TD', 'other_other'='OP_OP', 'other_PS'='OP_PS',
  'PS_TD'='PS_TD', 'PS_other'='PS_OP', 'PS_PS'='PS_PS')

final_df <- merge(clin_df, diag_df, by='bblid')

final_df$t1_tfinal <- relevel(final_df$t1_tfinal, 'TD_TD')

final_df$oT1_Tfinal <- ordered(final_df$t1_tfinal, c('TD_TD', 'OP_OP', 'OP_PS',
  'OP_TD', 'PS_OP', 'PS_PS', 'PS_TD', 'TD_OP', 'TD_PS'))

names(final_df)[names(final_df) == 'gaf_c'] <- 'GAF_Current'
names(final_df)[names(final_df) == 'gaf_h'] <- 'GAF_Highest'
names(final_df)[names(final_df) == 'age'] <- 'Age'

getUpperLowerCI <- function(i) {
  sorted_vec <- unname(sort(fits[i,]))
  lower <- sorted_vec[round(.025*length(sorted_vec))]
  upper <- sorted_vec[round(.975*length(sorted_vec))]
  c(lower, upper)
}

############################# Plot SIPS Dimensions #############################

plotcols <- c('Positive', 'Negative', 'Disorganized', 'General', 'GAF_Current', 'GAF_Highest')

diags <- as.character(unique(final_df$t1_tfinal)[!(unique(final_df$t1_tfinal) %in% c('TD_TD', 'PS_PS'))])

for (score in plotcols) {
  clin_score_df <- final_df
  names(clin_score_df)[names(clin_score_df) == score] <- 'score'
  clin_score_df <- clin_score_df[!is.na(clin_score_df$score), ]
  row.names(clin_score_df) <- 1:nrow(clin_score_df)

  #mod1_diag <- gamm4(score ~ t1_tfinal + s(Age, k=10, bs='cr'), data=clin_score_df,
  #  random=~(1|bblid), REML=TRUE)
  mod1b <- gamm4(score ~ t1_tfinal + s(Age, k=4, bs='cr') +
    s(Age, by=oT1_Tfinal, k=10, bs='cr'), data=clin_score_df, random=~(1|bblid), REML=TRUE)
  #tab_model(mod1_diag, mod2_diag)

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
    subtit <- paste0('Sessions: TD-TD=', nrow(clin_score_df[clin_score_df$t1_tfinal == 'TD_TD',]),
      ', ', group, '=', nrow(clin_score_df[clin_score_df$t1_tfinal == group,]),
      ', PS-PS=', nrow(clin_score_df[clin_score_df$t1_tfinal == 'PS_PS',]))
    group_df <- clin_score_df[clin_score_df$t1_tfinal %in% c('TD_TD', group, 'PS_PS'), ]
    group_df$t1_tfinal <- ordered(group_df$t1_tfinal, c('TD_TD', group, 'PS_PS'))
    clin_plot <- ggplot(group_df,
        aes(x=Age, y=score, color=t1_tfinal)) + theme_linedraw() +
      scale_color_manual(values=c('green3', 'goldenrod2', 'red')) +
      theme(legend.position = 'bottom', plot.title=element_text(size=14, face="bold"),
        plot.subtitle=element_text(size=8)) +
      labs(title=score, subtitle=subtit) + geom_line(aes(y=predgamm), size=1) +
      geom_line(data=clin_score_df[clin_score_df$t1_tfinal == 'TD_TD',], aes(y=LCI),
        size=.7, linetype=2, color='gray20') +
      geom_line(data=clin_score_df[clin_score_df$t1_tfinal == 'TD_TD',], aes(y=UCI),
        size=.7, linetype=2, color='gray20') +
      geom_line(data=clin_score_df[clin_score_df$t1_tfinal == group,], aes(y=LCI),
        size=.7, linetype=2, color='gray40') +
      geom_line(data=clin_score_df[clin_score_df$t1_tfinal == group,], aes(y=UCI),
        size=.7, linetype=2, color='gray40') +
      geom_line(data=clin_score_df[clin_score_df$t1_tfinal == 'PS_PS',], aes(y=LCI),
        size=.7, linetype=2, color='gray60') +
      geom_line(data=clin_score_df[clin_score_df$t1_tfinal == 'PS_PS',], aes(y=UCI),
        size=.7, linetype=2, color='gray60')

    assign(paste0(score, '_', group, '_plot'), clin_plot)

    pdf(file=paste0('~/Documents/pncLongitudinalPsychosis/plots/clin_', score, '_', group, '.pdf'), width=4, height=4)
    print(clin_plot)
    dev.off()
  }
}
