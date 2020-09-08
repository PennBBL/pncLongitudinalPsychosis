### This script does an exploratory factor analysis of the longitudinal CNB
### data, without accounting for repeated measures. In addition, it imputes
### missing data. It plots the data by sex, tests for sex differences within
### diagnostic groups, and includes MPRAXIS (speed) and TAP (speed) in the
### factor analysis.
###
### Ellyn Butler
### August 24, 2020 - August 27, 2020

set.seed(20)

library('dplyr')
library('reshape2')
library('ggplot2')
library('ggpubr')
library('psych')
library('missForest')
library('lme4')
library('gamm4')
library('stringr')
library('broom')
library('tidyr')
library('gratia')
library('MASS')

getUpperLowerCI <- function(i) {
  sorted_vec <- unname(sort(fits[i,]))
  lower <- sorted_vec[round(.025*length(sorted_vec))]
  upper <- sorted_vec[round(.975*length(sorted_vec))]
  c(lower, upper)
}

cnb_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/cnb_quickFA_impute_sex.csv')


###################### Plot factor scores ######################

plotcols <- c(paste0('ACC_Soln3_MR', 1:3), paste0('RT_Soln3_MR', 1:3),
  paste0('EFF_Soln4_MR', 1:4))

for (test in plotcols) {
  cnb_test_df <- cnb_df
  names(cnb_test_df)[names(cnb_test_df) == test] <- 'test'

  cnb_test_df$sex <- factor(cnb_test_df$sex)
  cnb_test_df$t1_tfinal_factor <- factor(cnb_test_df$t1_tfinal)

  #cnb_test_df$sex_t1_tfinal_factor <- paste0(cnb_test_df$sex, cnb_test_df$t1_tfinal_factor)

  mod1b <- gamm4(test ~ s(Age, k=10, bs='cr') + s(Age, by=sex, k=10, bs='cr') +
    s(Age, by=t1_tfinal_factor, k=10, bs='cr') +
    s(Age, by=sex:t1_tfinal_factor, k=10, bs='cr'),
    data=cnb_test_df, random=~(1|bblid), REML=TRUE)
    # August 25, 2020: Desired model matrix after dropping columns, but not sure
    # how to specify correctly so that dropping isn't needed. But t1_tfinal_factor
    # had to be an ordered factor for this to work
    # CONCLUSION: No interactions (so will run as separate models):
    # https://stats.stackexchange.com/questions/32730/how-to-include-an-interaction-term-in-gam
  #capture.output(gam.check(mod1b$gam),
  #  file=paste0('~/Documents/pncLongitudinalPsychosis/results/', test, '_sex_check_gamm_mod1b.txt'))

  for (lev in levels(cnb_test_df$t1_tfinal)) {
    dat <- cnb_test_df[cnb_test_df$t1_tfinal == lev, ]
    mod1b <- gamm4(test ~ s(Age, by=sex, k=10, bs='cr'), data=dat, random=~(1|bblid), REML=TRUE)
    # August 26, 2020: The main effect for age made this rank deficient...
    # but how do I interpret the p-values for the interaction then?

    lp <- predict(mod1b$gam, newdata=dat, type='lpmatrix')
    coefs <- coef(mod1b$gam) #Fits 262 coefficients
    vc <- vcov(mod1b$gam)

    sim <- mvrnorm(1000, mu = coefs, Sigma = vc)
    fits <- lp %*% t(sim)

    cis <- t(sapply(1:nrow(fits), getUpperLowerCI))
    cis <- data.frame(cis)
    names(cis) <- c('LCI', 'UCI')
    dat <- cbind(dat, cis)

    model_info <- tidy(mod1b$gam) %>%
      filter(str_detect(term, "sex"))
    #model_info <- model_info[model_info$p.value < .05,]
    model_info$t1_tfinal <- lev

    dat$predgamm <- predict(mod1b$gam)
    if (lev == levels(cnb_test_df$t1_tfinal)[1]) {
      cnb_test_df2 <- dat
    } else {
      cnb_test_df2 <- rbind(cnb_test_df2, dat)
    }
  }
  row.names(cnb_test_df2) <- 1:nrow(cnb_test_df2)
  cnb_test_df <- cnb_test_df2

  subtit <- paste0('Sessions: TD-TD=', nrow(cnb_test_df[cnb_test_df$t1_tfinal == 'TD_TD',]),
    ', TD-OP=', nrow(cnb_test_df[cnb_test_df$t1_tfinal == 'TD_OP',]),
    ', TD-PS=', nrow(cnb_test_df[cnb_test_df$t1_tfinal == 'TD_PS',]),
    ', OP-TD=', nrow(cnb_test_df[cnb_test_df$t1_tfinal == 'OP_TD',]),
    ', OP-OP=', nrow(cnb_test_df[cnb_test_df$t1_tfinal == 'OP_OP',]),
    ', OP-PS=', nrow(cnb_test_df[cnb_test_df$t1_tfinal == 'OP_PS',]),
    ', PS-TD=', nrow(cnb_test_df[cnb_test_df$t1_tfinal == 'PS_TD',]),
    ', PS-OP=', nrow(cnb_test_df[cnb_test_df$t1_tfinal == 'PS_OP',]),
    ', PS-PS=', nrow(cnb_test_df[cnb_test_df$t1_tfinal == 'PS_PS',]))

  cnb_test_df$first_diagnosis <- recode(cnb_test_df$first_diagnosis,
    'PS'='PS - First Diagnosis', 'OP'='OP - First Diagnosis', 'TD'='TD - First Diagnosis')
  cnb_test_df$first_diagnosis <- ordered(cnb_test_df$first_diagnosis,
    c('PS - First Diagnosis', 'OP - First Diagnosis', 'TD - First Diagnosis'))
  cnb_test_df$last_diagnosis <- recode(cnb_test_df$last_diagnosis,
    'TD'='TD - Last Diagnosis', 'OP'='OP - Last Diagnosis', 'PS'='PS - Last Diagnosis')
  cnb_test_df$last_diagnosis <- ordered(cnb_test_df$last_diagnosis,
    c('TD - Last Diagnosis', 'OP - Last Diagnosis', 'PS - Last Diagnosis'))

  cnb_plot <- ggplot(cnb_test_df, aes(x=Age, y=test, color=sex)) +
    theme_linedraw() + geom_line(aes(group=bblid), alpha=.2) +
    facet_grid(first_diagnosis ~ last_diagnosis) +
    scale_color_manual(values=c('red', 'blue')) +
    theme(legend.position = 'bottom', plot.title=element_text(size=14, face="bold"),
      plot.subtitle=element_text(size=8)) +
    labs(title=test, subtitle=subtit) + geom_line(aes(y=predgamm), size=1) +
    geom_line(data=cnb_test_df[cnb_test_df$sex == 'Female',], aes(y=LCI), size=.7,
      linetype=2, color='gray20') +
    geom_line(data=cnb_test_df[cnb_test_df$sex == 'Female',], aes(y=UCI), size=.7,
      linetype=2, color='gray20') +
    geom_line(data=cnb_test_df[cnb_test_df$sex == 'Male',], aes(y=LCI), size=.7,
      linetype=2, color='gray60') +
    geom_line(data=cnb_test_df[cnb_test_df$sex == 'Male',], aes(y=UCI), size=.7,
      linetype=2, color='gray60')

  #if (nrow(model_info) > 0) {
  #  # List the first-last pairs whose age trajectories significantly differ from TD-TD
  #  diagcats <- gsub('factor', '', model_info$term)
  #  diagcats <- gsub('s\\(Age\\):t1_tfinal_', '', diagcats)
  #  ann_text <- data.frame(t1_tfinal=diagcats, lab = "*", Age=28)
  #  int <- strsplit(as.character(ann_text$t1_tfinal), '_')
  #  ann_text$first_diagnosis <- sapply(int, `[[`, 1)
  #  ann_text$first_diagnosis <- paste(ann_text$first_diagnosis, '- First Diagnosis')
  #  ann_text$last_diagnosis <- sapply(int, `[[`, 2)
  #  ann_text$last_diagnosis <- paste(ann_text$last_diagnosis, '- Last Diagnosis')

  #  ann_text$first_diagnosis <- ordered(ann_text$first_diagnosis,
  #    c('PS - First Diagnosis', 'OP - First Diagnosis', 'TD - First Diagnosis'))
  #  ann_text$last_diagnosis <- recode(ann_text$last_diagnosis,
  #    'TD'='TD - Last Diagnosis', 'OP'='OP - Last Diagnosis', 'PS'='PS - Last Diagnosis')
  #  ann_text$last_diagnosis <- ordered(ann_text$last_diagnosis,
  #    c('TD - Last Diagnosis', 'OP - Last Diagnosis', 'PS - Last Diagnosis'))

  #  num <- min(cnb_test_df$test) + (max(cnb_test_df$test) - min(cnb_test_df$test))*.33
  #  cnb_plot <- cnb_plot + geom_text(data=ann_text, y=num, label="*", size=15)
  #}
  assign(paste0(test, '_plot'), cnb_plot)
}


pdf(file='~/Documents/pncLongitudinalPsychosis/plots/longFactor3x3_impute_sex.pdf', width=7.5, height=6)
ACC_Soln3_MR1_plot
ACC_Soln3_MR2_plot
ACC_Soln3_MR3_plot
RT_Soln3_MR1_plot
RT_Soln3_MR2_plot
RT_Soln3_MR3_plot
EFF_Soln4_MR1_plot
EFF_Soln4_MR2_plot
EFF_Soln4_MR3_plot
EFF_Soln4_MR4_plot
dev.off()
