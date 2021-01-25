### This script plots the fits for WIV (accuracy)
###
### Ellyn Butler
### December 16, 2020

set.seed(20)

library('dplyr')
library('reshape2')
library('ggplot2')
library('ggpubr')
library('gamm4')
library('stringr')
library('MASS')
library('broom')
library('sjPlot')

final_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/cognitive/cnb_quickFA_impute_2020-10-20.csv', stringsAsFactors=TRUE)

final_df[, grep('_ACC', names(final_df), value=TRUE)] <- sapply(final_df[, grep('_ACC', names(final_df), value=TRUE)], scale)
final_df$WIV <- scale(apply(final_df[, grep('_ACC', names(final_df), value=TRUE)], 1, sd))

final_df$t1_tfinal <- recode(final_df$t1_tfinal, 'TD_TD'='TD-TD', 'OP_OP'='OP-OP',
  'OP_PS'='OP-PS', 'OP_TD'='OP-TD', 'PS_OP'='PS-OP', 'PS_PS'='PS-PS', 'PS_TD'='PS-TD',
  'TD_OP'='TD-OP', 'TD_PS'='TD-PS')
final_df$t1_tfinal <- relevel(final_df$t1_tfinal, 'TD-TD')

final_df$oT1_Tfinal <- ordered(final_df$t1_tfinal, c('TD-TD', 'OP-OP', 'OP-PS',
  'OP-TD', 'PS-OP', 'PS-PS', 'PS-TD', 'TD-OP', 'TD-PS'))


getUpperLowerCI <- function(i) {
  sorted_vec <- unname(sort(fits[i,]))
  lower <- sorted_vec[round(.025*length(sorted_vec))]
  upper <- sorted_vec[round(.975*length(sorted_vec))]
  c(lower, upper)
}

###################### Plot WIV ######################


diags <- as.character(unique(final_df$t1_tfinal)[!(unique(final_df$t1_tfinal) %in% c('TD-TD', 'PS-PS'))])

final_df$t1_tfinal_factor <- factor(final_df$t1_tfinal)

mod1b <- gamm4(WIV ~ t1_tfinal +  s(Age, k=4, bs='cr') + s(Age, by=oT1_Tfinal, k=4, bs='cr'),
  data=final_df, random=~(1|bblid), REML=TRUE)

lp <- predict(mod1b$gam, newdata=final_df, type='lpmatrix')
coefs <- coef(mod1b$gam)
vc <- vcov(mod1b$gam)

sim <- mvrnorm(1000, mu = coefs, Sigma = vc)
fits <- lp %*% t(sim)

cis <- t(sapply(1:nrow(fits), getUpperLowerCI))
cis <- data.frame(cis)
names(cis) <- c('LCI', 'UCI')
final_df <- cbind(final_df, cis)

model_info <- tidy(mod1b$gam) %>%
  filter(str_detect(term, "t1_tfinal_factor"))

model_info <- model_info[model_info$p.value < .05,]

final_df$predgamm <- predict(mod1b$gam)

#tab_model(mod1b$gam)

for (group in diags) {
    subtit <- paste0('Sessions: TD-TD=', nrow(final_df[final_df$t1_tfinal == 'TD-TD',]),
      ', ', group, '=', nrow(final_df[final_df$t1_tfinal == group,]),
      ', PS-PS=', nrow(final_df[final_df$t1_tfinal == 'PS-PS',]))
    group_df <- final_df[final_df$t1_tfinal %in% c('TD-TD', group, 'PS-PS'), ]
    group_df$t1_tfinal <- ordered(group_df$t1_tfinal, c('TD-TD', group, 'PS-PS'))
    cnb_plot <- ggplot(group_df,
        aes(x=Age, y=WIV, color=t1_tfinal)) + theme_linedraw() +
      scale_color_manual(values=c('green3', 'goldenrod2', 'red')) +
      theme(legend.position = 'bottom', plot.title=element_text(size=14, face="bold"),
        plot.subtitle=element_text(size=8)) +
      labs(title='Accuracy WIV', subtitle=subtit) + geom_line(aes(y=predgamm), size=1) +
      geom_line(data=final_df[final_df$t1_tfinal == 'TD-TD',], aes(y=LCI),
        size=.7, linetype=2, color='gray20') +
      geom_line(data=final_df[final_df$t1_tfinal == 'TD-TD',], aes(y=UCI),
        size=.7, linetype=2, color='gray20') +
      geom_line(data=final_df[final_df$t1_tfinal == group,], aes(y=LCI),
        size=.7, linetype=2, color='gray40') +
      geom_line(data=final_df[final_df$t1_tfinal == group,], aes(y=UCI),
        size=.7, linetype=2, color='gray40') +
      geom_line(data=final_df[final_df$t1_tfinal == 'PS-PS',], aes(y=LCI),
        size=.7, linetype=2, color='gray60') +
      geom_line(data=final_df[final_df$t1_tfinal == 'PS-PS',], aes(y=UCI),
        size=.7, linetype=2, color='gray60')

    pdf(file=paste0('~/Documents/pncLongitudinalPsychosis/plots/cnbFactorImpute_ACC_WIV_', group, '.pdf'), width=4, height=4)
    print(cnb_plot)
    dev.off()
}












#
