### This script plots just the fits for the factors derived in quickFAImpute.R
###
### Ellyn Butler
### August 27, 2020 - September 21, 2020


library('dplyr')
library('reshape2')
library('ggplot2')
library('ggpubr')
library('lme4')
library('gamm4')
library('stringr')
library('MASS')
library('broom')

final_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/cognitive/cnb_quickFA_impute_2020-10-20.csv')

final_df$t1_tfinal <- ordered(final_df$t1_tfinal, c('TD_TD', 'TD_OP', 'TD_PS',
  'OP_TD', 'OP_OP', 'OP_PS', 'PS_TD', 'PS_OP', 'PS_PS'))

names(final_df)[names(final_df) == 'EFF_Soln4_MR1'] <- 'SocCog_EFF'
names(final_df)[names(final_df) == 'EFF_Soln4_MR2'] <- 'Exec_EFF'
names(final_df)[names(final_df) == 'EFF_Soln4_MR3'] <- 'Mem_EFF'
names(final_df)[names(final_df) == 'EFF_Soln4_MR4'] <- 'CompCog_EFF'

final_df$oSex <- ordered(final_df$sex, c('Male', 'Female'))

getUpperLowerCI <- function(i) {
  sorted_vec <- unname(sort(fits[i,]))
  lower <- sorted_vec[round(.025*length(sorted_vec))]
  upper <- sorted_vec[round(.975*length(sorted_vec))]
  c(lower, upper)
}

###################### Plot factor scores ######################

plotcols <- c('SocCog_EFF', 'Exec_EFF', 'Mem_EFF', 'CompCog_EFF')

for (test in plotcols) {
  cnb_test_df <- final_df
  names(cnb_test_df)[names(cnb_test_df) == test] <- 'test'

  mod1b <- gamm4(test ~ sex +  s(Age, by=oSex, k=10, bs='cr') +
    s(Age, k=10, bs='cr'), data=cnb_test_df, random=~(1|bblid), REML=TRUE)

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

  subtit <- paste0('Sessions: Female=', nrow(cnb_test_df[cnb_test_df$sex == 'Female',]),
    ', Male=', nrow(cnb_test_df[cnb_test_df$sex == 'Male',]))

  cnb_plot <- ggplot(cnb_test_df, aes(x=Age, y=test, color=sex)) + theme_linedraw() +
    scale_color_manual(values=c('red', 'blue')) +
    theme(legend.position = 'bottom', plot.title=element_text(size=14, face="bold"),
      plot.subtitle=element_text(size=8)) +
    labs(title=test, subtitle=subtit) + geom_line(aes(y=predgamm), size=1) +
    geom_line(data=cnb_test_df[cnb_test_df$sex == 'Female',], aes(y=LCI),
      size=.7, linetype=2, color='gray20') +
    geom_line(data=cnb_test_df[cnb_test_df$sex == 'Female',], aes(y=UCI),
      size=.7, linetype=2, color='gray20') +
    geom_line(data=cnb_test_df[cnb_test_df$sex == 'Male',], aes(y=LCI),
      size=.7, linetype=2, color='gray60') +
    geom_line(data=cnb_test_df[cnb_test_df$sex== 'Male',], aes(y=UCI),
      size=.7, linetype=2, color='gray60')

  assign(paste0(test, '_plot'), cnb_plot)
}

pdf(file='~/Documents/pncLongitudinalPsychosis/plots/cnbFactorImpute_SocCog_EFF_sex.pdf', width=4, height=4)
SocCog_EFF_plot
dev.off()

pdf(file='~/Documents/pncLongitudinalPsychosis/plots/cnbFactorImpute_Exec_EFF_sex.pdf', width=4, height=4)
Exec_EFF_plot
dev.off()

pdf(file='~/Documents/pncLongitudinalPsychosis/plots/cnbFactorImpute_Mem_EFF_sex.pdf', width=4, height=4)
Mem_EFF_plot
dev.off()

pdf(file='~/Documents/pncLongitudinalPsychosis/plots/cnbFactorImpute_CompCog_EFF_sex.pdf', width=4, height=4)
CompCog_EFF_plot
dev.off()
