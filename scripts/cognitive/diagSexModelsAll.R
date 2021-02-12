### This script runs GAMMs to see if there are main effects and interactions
### for sex and longitudinal diagnostic labels for all of the seven WIV variables
###
### Ellyn Butler
### January 11, 2021 - January 27, 2021


set.seed(20)

library('dplyr')
library('reshape2')
library('ggpubr')
library('psych')
library('lme4')
library('gamm4')
library('sjPlot')
library('MASS')

cnb_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/cognitive/cnb_quickFA_impute_2020-10-20.csv', stringsAsFactors = TRUE)

# Recode race
cnb_df$race <- recode(cnb_df$race, `1`='White', `2`='Other', `3`='Other',
  `4`='Other', `5`='Other')

############## Create the WIV metrics
cnb_df[, grep('_EFF', names(cnb_df), value=TRUE)] <- sapply(cnb_df[, grep('_EFF', names(cnb_df), value=TRUE)], scale)
cnb_df$WIV_EFF <- scale(apply(cnb_df[, grep('_EFF', names(cnb_df), value=TRUE)], 1, sd))

cnb_df[, grep('_ACC', names(cnb_df), value=TRUE)] <- sapply(cnb_df[, grep('_ACC', names(cnb_df), value=TRUE)], scale)
cnb_df$WIV_ACC <- scale(apply(cnb_df[, grep('_ACC', names(cnb_df), value=TRUE)], 1, sd))

cnb_df[, grep('_RT', names(cnb_df), value=TRUE)] <- sapply(cnb_df[, grep('_RT', names(cnb_df), value=TRUE)], scale)
cnb_df$WIV_RT <- scale(apply(cnb_df[, grep('_RT', names(cnb_df), value=TRUE)], 1, sd))

cnb_df$WIV_Exec_EFF <- scale(apply(cnb_df[, c('CPT_EFF', 'NBACK_EFF', 'PCET_EFF')], 1, sd))
cnb_df$WIV_Mem_EFF <- scale(apply(cnb_df[, c('CPF_EFF', 'CPW_EFF', 'VOLT_EFF')], 1, sd))
cnb_df$WIV_CompCog_EFF <- scale(apply(cnb_df[, c('PLOT_EFF', 'PMAT_EFF', 'PVRT_EFF')], 1, sd))
cnb_df$WIV_SocCog_EFF <- scale(apply(cnb_df[, c('ADT_EFF', 'ER40_EFF', 'MEDF_EFF')], 1, sd))
# ^ NOTE: Included PCET with Exec, even though it loads much higher on the CompCog
# factor, because historically it has loaded on the Exec factor and we need at
# least three tests to compute a SD

############## Name and scale factor scores
names(cnb_df)[names(cnb_df) == 'EFF_Soln4_MR1'] <- 'SocCog_EFF'
names(cnb_df)[names(cnb_df) == 'EFF_Soln4_MR2'] <- 'Exec_EFF'
names(cnb_df)[names(cnb_df) == 'EFF_Soln4_MR3'] <- 'Mem_EFF'
names(cnb_df)[names(cnb_df) == 'EFF_Soln4_MR4'] <- 'CompCog_EFF'

cnb_df$SocCog_EFF <- scale(cnb_df$SocCog_EFF)
cnb_df$Exec_EFF <- scale(cnb_df$Exec_EFF)
cnb_df$Mem_EFF <- scale(cnb_df$Mem_EFF)
cnb_df$CompCog_EFF <- scale(cnb_df$CompCog_EFF)

############## Recoding variables
cnb_df$sex <- relevel(cnb_df$sex, 'Male')
cnb_df$t1_tfinal <- relevel(cnb_df$t1_tfinal, 'TD_TD')

cnb_df$oSex <- ordered(cnb_df$sex, c('Male', 'Female'))

cnb_df$t1_tfinal <- recode(cnb_df$t1_tfinal, 'TD_TD'='TD-TD', 'OP_OP'='OP-OP',
  'OP_PS'='OP-PS', 'OP_TD'='OP-TD', 'PS_OP'='PS-OP', 'PS_PS'='PS-PS', 'PS_TD'='PS-TD',
  'TD_OP'='TD-OP', 'TD_PS'='TD-PS')

cnb_df$oT1_Tfinal <- ordered(cnb_df$t1_tfinal, c('TD-TD', 'OP-OP', 'OP-PS',
  'OP-TD', 'PS-OP', 'PS-PS', 'PS-TD', 'TD-OP', 'TD-PS'))

cnb_df$sex_t1_tfinal <- paste(cnb_df$sex, cnb_df$t1_tfinal, sep='_')
cnb_df$oSex_oT1_Tfinal <- ordered(cnb_df$sex_t1_tfinal, c('Male_TD-TD',
  'Male_OP-OP', 'Male_OP-PS', 'Male_OP-TD', 'Male_PS-OP', 'Male_PS-PS',
  'Male_PS-TD', 'Male_TD-OP', 'Male_TD-PS', 'Female_TD-TD', 'Female_OP-OP',
  'Female_OP-PS', 'Female_OP-TD', 'Female_PS-OP', 'Female_PS-PS', 'Female_PS-TD',
  'Female_TD-OP', 'Female_TD-PS'))

############## Defining functions
getUpperLowerCI <- function(i) {
  sorted_vec <- unname(sort(fits[i,]))
  lower <- sorted_vec[round(.025*length(sorted_vec))]
  upper <- sorted_vec[round(.975*length(sorted_vec))]
  c(lower, upper)
}

################################## Titles ##################################

varis <- c(grep('WIV', names(cnb_df), value=TRUE), 'SocCog_EFF', 'Exec_EFF', 'Mem_EFF', 'CompCog_EFF')

title_df <- data.frame(Variable=varis, NewName=c('WIV Efficiency', 'WIV Accuracy',
  'WIV Speed', 'WIV Executive Efficiency', 'WIV Memory Efficiency', 'WIV Complex Efficiency',
  'WIV Social Efficiency', 'Social Efficiency', 'Executive Efficiency',
  'Memory Efficiency', 'Complex Efficiency'), LB=c(rep(-1.5, 3), rep(-1.5, 4),
  rep(-3, 4)), UB=c(rep(2, 11)))


################################## GAMMs ##################################

for (vari in varis) {
  mod1 <- gamm4(as.formula(paste(vari , "~ t1_tfinal + s(Age, k=4, bs='cr') +
    s(Age, by=oT1_Tfinal, k=4, bs='cr')")), data=cnb_df, random=~(1|bblid), REML=TRUE)
  mod2 <- gamm4(as.formula(paste(vari , "~ sex + race + t1_tfinal + s(Age, k=4, bs='cr') +
    s(Age, by=oT1_Tfinal, k=4, bs='cr')")), data=cnb_df, random=~(1|bblid), REML=TRUE)
  print(tab_model(mod1$gam, mod2$gam, file=paste0('~/Documents/pncLongitudinalPsychosis/results/table_', vari, '.html')))

  # Get plot info
  LBY <- title_df[title_df$Variable == vari, 'LB']
  UBY <- title_df[title_df$Variable == vari, 'UB']
  TIT <- title_df[title_df$Variable == vari, 'NewName']

  ### Split by sex for plotting, and not
  for (sex in c(as.character(unique(cnb_df$sex)), 'both')) {
    if (sex == 'both') { final_df <- cnb_df } else { final_df <- cnb_df[cnb_df$sex == sex,] }
    row.names(final_df) <- 1:nrow(final_df)
    diags <- as.character(unique(final_df$t1_tfinal)[!(unique(final_df$t1_tfinal) %in% c('TD-TD', 'PS-PS'))])

    final_df$t1_tfinal_factor <- factor(final_df$t1_tfinal)

    mod1b <- gamm4(as.formula(paste(vari , "~ t1_tfinal +  s(Age, k=4, bs='cr') +
      s(Age, by=oT1_Tfinal, k=4, bs='cr')")), data=final_df, random=~(1|bblid), REML=TRUE)

    print(tab_model(mod1b$gam, file=paste0('~/Documents/pncLongitudinalPsychosis/results/table_', vari, '_', sex, '.html'))) # >>> Model Sex Split

    lp <- predict(mod1b$gam, newdata=final_df, type='lpmatrix')
    coefs <- coef(mod1b$gam)
    vc <- vcov(mod1b$gam)

    sim <- mvrnorm(1000, mu = coefs, Sigma = vc)
    fits <- lp %*% t(sim)

    cis <- t(sapply(1:nrow(fits), getUpperLowerCI))
    cis <- data.frame(cis)
    names(cis) <- c('LCI', 'UCI')
    final_df <- cbind(final_df, cis)
    final_df$predgamm <- predict(mod1b$gam)

    for (group in diags) {
        subtit <- paste0('Sessions: TD-TD=', nrow(final_df[final_df$t1_tfinal == 'TD-TD',]),
          ', ', group, '=', nrow(final_df[final_df$t1_tfinal == group,]),
          ', PS-PS=', nrow(final_df[final_df$t1_tfinal == 'PS-PS',]))
        group_df <- final_df[final_df$t1_tfinal %in% c('TD-TD', group, 'PS-PS'), ]
        group_df$t1_tfinal <- ordered(group_df$t1_tfinal, c('TD-TD', group, 'PS-PS'))
        cnb_plot <- ggplot(group_df, aes_string(x='Age', y=vari, color='t1_tfinal')) + theme_linedraw() + #Plot Base and Sex Split
          scale_color_manual(values=c('green3', 'goldenrod2', 'red')) +
          theme(legend.position = 'bottom', plot.title=element_text(size=14, face="bold"),
            plot.subtitle=element_text(size=8)) +
          ylim(LBY, UBY) + xlim(5, 30) + geom_hline(yintercept=0, size=1.5) + 
          labs(title=paste0(TIT, ' (', sex, ')'), subtitle=subtit, y='Score (95% CI)', color='Diagnoses') +
          geom_line(aes(y=predgamm), size=1) +
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

        pdf(file=paste0('~/Documents/pncLongitudinalPsychosis/plots/cnbFactorImpute_', vari, '_', group, '_', sex, '.pdf'), width=4, height=4)
        print(cnb_plot)
        dev.off()
    }
  }
}
