### This script plots just the fits for the factors derived in quickFAImpute.R
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

final_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/cognitive/cnb_quickFA_impute_2020-10-20.csv')

final_df$t1_tfinal <- relevel(final_df$t1_tfinal, 'TD_TD')

names(final_df)[names(final_df) == 'EFF_Soln4_MR1'] <- 'SocCog_EFF'
names(final_df)[names(final_df) == 'EFF_Soln4_MR2'] <- 'Exec_EFF'
names(final_df)[names(final_df) == 'EFF_Soln4_MR3'] <- 'Mem_EFF'
names(final_df)[names(final_df) == 'EFF_Soln4_MR4'] <- 'CompCog_EFF'

final_df$oT1_Tfinal <- ordered(final_df$t1_tfinal, c('TD_TD', 'OP_OP', 'OP_PS',
  'OP_TD', 'PS_OP', 'PS_PS', 'PS_TD', 'TD_OP', 'TD_PS'))


getUpperLowerCI <- function(i) {
  sorted_vec <- unname(sort(fits[i,]))
  lower <- sorted_vec[round(.025*length(sorted_vec))]
  upper <- sorted_vec[round(.975*length(sorted_vec))]
  c(lower, upper)
}

###################### Plot factor scores ######################

plotcols <- c('SocCog_EFF', 'Exec_EFF', 'Mem_EFF', 'CompCog_EFF')

diags <- as.character(unique(final_df$t1_tfinal)[!(unique(final_df$t1_tfinal) %in% c('TD_TD', 'PS_PS'))])

for (test in plotcols) {
  cnb_test_df <- final_df
  names(cnb_test_df)[names(cnb_test_df) == test] <- 'test'

  cnb_test_df$t1_tfinal_factor <- factor(cnb_test_df$t1_tfinal)

  mod1b <- gamm4(test ~ t1_tfinal +  s(Age, by=oT1_Tfinal, k=4, bs='cr') +
    s(Age, k=4, bs='cr'), data=cnb_test_df, random=~(1|bblid), REML=TRUE)
    # August 27, 2020: Ordering t1_tfinal fixes the need to drop a coefficient

  lp <- predict(mod1b$gam, newdata=cnb_test_df, type='lpmatrix')
  coefs <- coef(mod1b$gam)
  vc <- vcov(mod1b$gam)

  sim <- mvrnorm(1000, mu = coefs, Sigma = vc)
  fits <- lp %*% t(sim)

  cis <- t(sapply(1:nrow(fits), getUpperLowerCI))
  cis <- data.frame(cis)
  names(cis) <- c('LCI', 'UCI')
  cnb_test_df <- cbind(cnb_test_df, cis)

  model_info <- tidy(mod1b$gam) %>%
    filter(str_detect(term, "t1_tfinal_factor"))
  assign(paste0(test, '_model'), model_info)
  model_info <- model_info[model_info$p.value < .05,]

  cnb_test_df$predgamm <- predict(mod1b$gam)

  for (group in diags) {
    subtit <- paste0('Sessions: TD-TD=', nrow(cnb_test_df[cnb_test_df$t1_tfinal == 'TD_TD',]),
      ', ', group, '=', nrow(cnb_test_df[cnb_test_df$t1_tfinal == group,]),
      ', PS-PS=', nrow(cnb_test_df[cnb_test_df$t1_tfinal == 'PS_PS',]))
    group_df <- cnb_test_df[cnb_test_df$t1_tfinal %in% c('TD_TD', group, 'PS_PS'), ]
    group_df$t1_tfinal <- ordered(group_df$t1_tfinal, c('TD_TD', group, 'PS_PS'))
    cnb_plot <- ggplot(group_df,
        aes(x=Age, y=test, color=t1_tfinal)) + theme_linedraw() +
      scale_color_manual(values=c('green3', 'goldenrod2', 'red')) +
      theme(legend.position = 'bottom', plot.title=element_text(size=14, face="bold"),
        plot.subtitle=element_text(size=8)) +
      labs(title=test, subtitle=subtit) + geom_line(aes(y=predgamm), size=1) +
      geom_line(data=cnb_test_df[cnb_test_df$t1_tfinal == 'TD_TD',], aes(y=LCI),
        size=.7, linetype=2, color='gray20') +
      geom_line(data=cnb_test_df[cnb_test_df$t1_tfinal == 'TD_TD',], aes(y=UCI),
        size=.7, linetype=2, color='gray20') +
      geom_line(data=cnb_test_df[cnb_test_df$t1_tfinal == group,], aes(y=LCI),
        size=.7, linetype=2, color='gray40') +
      geom_line(data=cnb_test_df[cnb_test_df$t1_tfinal == group,], aes(y=UCI),
        size=.7, linetype=2, color='gray40') +
      geom_line(data=cnb_test_df[cnb_test_df$t1_tfinal == 'PS_PS',], aes(y=LCI),
        size=.7, linetype=2, color='gray60') +
      geom_line(data=cnb_test_df[cnb_test_df$t1_tfinal == 'PS_PS',], aes(y=UCI),
        size=.7, linetype=2, color='gray60')

    assign(paste0(test, '_', group, '_plot'), cnb_plot)

    pdf(file=paste0('~/Documents/pncLongitudinalPsychosis/plots/cnbFactorImpute_', test, '_', group, '.pdf'), width=4, height=4)
    print(cnb_plot)
    dev.off()
  }
}











#
