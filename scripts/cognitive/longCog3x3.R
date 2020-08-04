### This script plots longitudinal trajectories of CNB tests by first and last
### diagnosis
###
### Ellyn Butler
### July 27, 2020 - August 4, 2020

set.seed(20)

library('dplyr')
library('reshape2')
library('ggplot2')
library('ggpubr')
library('lme4')
library('gamm4')
library('stringr')
library('broom')
library('tidyr')

clin_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/clinical/pnc_longitudinal_diagnosis_n752_202007.csv')
names(clin_df)[names(clin_df) == 't1'] <- 'first_diagnosis'
names(clin_df)[names(clin_df) == 'tfinal2'] <- 'last_diagnosis'
clin_df$first_diagnosis <- as.character(clin_df$first_diagnosis)
clin_df$last_diagnosis <- as.character(clin_df$last_diagnosis)
clin_df$first_diagnosis <- recode(clin_df$first_diagnosis, 'TD'='TD', 'other'='OP', 'PS'='PS')
clin_df$last_diagnosis <- recode(clin_df$last_diagnosis, 'TD'='TD', 'other'='OP', 'PS'='PS')
clin_df$t1_tfinal <- recode(clin_df$t1_tfinal, 'TD_TD'='TD_TD', 'TD_other'='TD_OP',
  'TD_PS'='TD_PS', 'other_TD'='OP_TD', 'other_other'='OP_OP', 'other_PS'='OP_PS',
  'PS_TD'='PS_TD', 'PS_other'='PS_OP', 'PS_PS'='PS_PS')

cnb_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/cognitive/CNB_Longitudinal_Core_11February2020.csv')
cnb_df <- cnb_df[, c('bblid', 'Age', 'timepoint', 'Test', 'ACC_raw')]
cnb_df <- cnb_df[cnb_df$bblid %in% clin_df$bblid,]
cnb_df <- dcast(cnb_df, bblid + Age + timepoint ~ Test, value.var='ACC_raw')

getDiagnoses <- function(i) {
  bblid <- cnb_df[i, 'bblid']
  c(clin_df[clin_df$bblid == bblid, 'first_diagnosis'], clin_df[clin_df$bblid == bblid, 'last_diagnosis'], as.character(clin_df[clin_df$bblid == bblid, 't1_tfinal']))
}

cnb_df[,c('first_diagnosis', 'last_diagnosis', 't1_tfinal')] <- t(sapply(1:nrow(cnb_df), getDiagnoses))

cnb_df$first_diagnosis <- ordered(cnb_df$first_diagnosis, c('TD', 'OP', 'PS'))
cnb_df$last_diagnosis <- ordered(cnb_df$last_diagnosis, c('TD', 'OP', 'PS'))
cnb_df$t1_tfinal <- ordered(cnb_df$t1_tfinal, c('TD_TD', 'TD_OP', 'TD_PS',
  'OP_TD', 'OP_OP', 'OP_PS', 'PS_TD', 'PS_OP', 'PS_PS'))

for (test in c('ADT', 'CPF', 'CPT', 'CPW', 'ER40', 'MEDF', 'NBACK', 'PCET',
  'PLOT', 'PMAT', 'PVRT', 'TAP', 'VOLT')) {
  cnb_test_df <- cnb_df[!is.na(cnb_df[,test]),]
  row.names(cnb_test_df) <- 1:nrow(cnb_test_df)
  names(cnb_test_df)[names(cnb_test_df) == test] <- 'test'

  cnb_test_df$t1_tfinal_factor <- factor(cnb_test_df$t1_tfinal)

  mod1b <- gamm4(test ~ s(Age, by=t1_tfinal_factor, k=60, bs="cr") +
    s(Age, k=40, bs="cr"), data=cnb_test_df, random=~(1|bblid), REML=TRUE)
    # August 4: Without the second term, just fitting an intercept for TD-TD
  capture.output(gam.check(mod1b$gam),
    file=paste0('~/Documents/pncLongitudinalPsychosis/results/', test, '_check_gamm_mod1b.txt'))
  # ^ Can't get k-index above 1 (up to k=50), but k' is very far away from edf
  #TO DO: Save check plots

  #mod2b <- gamm4(test ~ t2(Age_bl, Time, k=c(20, 5), bs='cr'), data=cnb_test_df,
  #  random=~(1|bblid)) # Need to calculate Age_bl and Time, and then figure out the fits

  model_info <- tidy(mod1b$gam) %>%
    filter(str_detect(term, "t1_tfinal_factor"))
  assign(paste0(test, '_model'), model_info)
  model_info <- model_info[model_info$p.value < .05,]

  cnb_test_df$predgamm <- predict(mod1b$gam)

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

  cnb_plot <- ggplot(cnb_test_df, aes(x=Age, y=test, color=last_diagnosis)) +
    theme_linedraw() + geom_line(aes(group=bblid), alpha=.2) +
    facet_grid(first_diagnosis ~ last_diagnosis) +
    scale_color_manual(values=c('#009E73', '#CC79A7', '#0072B2')) +
    theme(legend.position = 'none', plot.title=element_text(size=14, face="bold"),
      plot.subtitle=element_text(size=8)) +
    labs(title=test, subtitle=subtit) + geom_line(aes(y=predgamm), size=1)

  if (nrow(model_info) > 0) {
    # List the first-last pairs whose age trajectories significantly differ from TD-TD
    diagcats <- gsub('factor', '', model_info$term)
    diagcats <- gsub('s\\(Age\\):t1_tfinal_', '', diagcats)
    ann_text <- data.frame(t1_tfinal=diagcats, lab = "*", Age=28)
    int <- strsplit(as.character(ann_text$t1_tfinal), '_')
    ann_text$first_diagnosis <- sapply(int, `[[`, 1)
    ann_text$first_diagnosis <- paste(ann_text$first_diagnosis, '- First Diagnosis')
    ann_text$last_diagnosis <- sapply(int, `[[`, 2)
    ann_text$last_diagnosis <- paste(ann_text$last_diagnosis, '- Last Diagnosis')

    ann_text$first_diagnosis <- ordered(ann_text$first_diagnosis,
      c('PS - First Diagnosis', 'OP - First Diagnosis', 'TD - First Diagnosis'))
    ann_text$last_diagnosis <- recode(ann_text$last_diagnosis,
      'TD'='TD - Last Diagnosis', 'OP'='OP - Last Diagnosis', 'PS'='PS - Last Diagnosis')
    ann_text$last_diagnosis <- ordered(ann_text$last_diagnosis,
      c('TD - Last Diagnosis', 'OP - Last Diagnosis', 'PS - Last Diagnosis'))

    num <- min(cnb_test_df$test) + (max(cnb_test_df$test) - min(cnb_test_df$test))*.33
    cnb_plot <- cnb_plot + geom_text(data=ann_text, y=num, label="*", size=15)
  }
  assign(paste0(test, '_plot'), cnb_plot)
}

# TO DO:
# https://stackoverflow.com/questions/31075407/plot-mixed-effects-model-in-ggplot



pdf(file='~/Documents/pncLongitudinalPsychosis/plots/longCog3x3.pdf', width=7.5, height=6)
ADT_plot
CPF_plot
CPT_plot
CPW_plot
ER40_plot
MEDF_plot
NBACK_plot
PCET_plot
PLOT_plot
PMAT_plot
PVRT_plot
TAP_plot
VOLT_plot
dev.off()










# Baseline groups
#model_first <- gamm4(test ~ s(Age, k=20, bs="cr") + s(Age, by=first_diagnosis, k=10, bs="cr"),
#  data=cnb_test_df, random=~(1|bblid))
#tidy(model_first$gam) %>%
#  filter(str_detect(term, "first_diagnosis"))

# Baseline groups
#model_last <- gamm4(test ~ s(Age, k=20, bs="cr") + s(Age, by=last_diagnosis, k=10, bs="cr"),
#  data=cnb_test_df, random=~(1|bblid))
#tidy(model_last$gam) %>%
#  filter(str_detect(term, "last_diagnosis"))
