### This script plots longitudinal trajectories of CNB tests by first and last
### diagnosis
###
### Ellyn Butler
### July 27, 2020


library('dplyr')
library('reshape2')
library('ggplot2')
library('ggpubr')

clin_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/clinical/pnc_longitudinal_diagnosis_n752_202007.csv')
names(clin_df)[names(clin_df) == 't1'] <- 'first_diagnosis'
names(clin_df)[names(clin_df) == 'tfinal2'] <- 'last_diagnosis'
clin_df$first_diagnosis <- as.character(clin_df$first_diagnosis)
clin_df$last_diagnosis <- as.character(clin_df$last_diagnosis)
clin_df$first_diagnosis <- recode(clin_df$first_diagnosis, 'TD'='TD', 'other'='OP', 'PS'='PS')
clin_df$last_diagnosis <- recode(clin_df$last_diagnosis, 'TD'='TD', 'other'='OP', 'PS'='PS')

cnb_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/cognitive/CNB_Longitudinal_Core_11February2020.csv')
cnb_df <- cnb_df[, c('bblid', 'Age', 'timepoint', 'Test', 'ACC_raw')]
cnb_df <- cnb_df[cnb_df$bblid %in% clin_df$bblid,]
cnb_df <- dcast(cnb_df, bblid + Age + timepoint ~ Test, value.var='ACC_raw')

getDiagnoses <- function(i) {
  bblid <- cnb_df[i, 'bblid']
  c(clin_df[clin_df$bblid == bblid, 'first_diagnosis'], clin_df[clin_df$bblid == bblid, 'last_diagnosis'])
}

cnb_df[,c('first_diagnosis', 'last_diagnosis')] <- t(sapply(1:nrow(cnb_df), getDiagnoses))

cnb_df$first_diagnosis <- recode(cnb_df$first_diagnosis,
  'PS'='PS - First Diagnosis', 'OP'='OP - First Diagnosis', 'TD'='TD - First Diagnosis')
cnb_df$first_diagnosis <- ordered(cnb_df$first_diagnosis,
  c('PS - First Diagnosis', 'OP - First Diagnosis', 'TD - First Diagnosis'))
cnb_df$last_diagnosis <- recode(cnb_df$last_diagnosis,
  'TD'='TD - Last Diagnosis', 'OP'='OP - Last Diagnosis', 'PS'='PS - Last Diagnosis')
cnb_df$last_diagnosis <- ordered(cnb_df$last_diagnosis,
  c('TD - Last Diagnosis', 'OP - Last Diagnosis', 'PS - Last Diagnosis'))


for (test in c('ADT', 'CPF', 'CPT', 'CPW', 'ER40', 'MEDF', 'MPRAXIS', 'NBACK', 'PCET',
  'PLOT', 'PMAT', 'PVRT', 'TAP', 'VOLT')) {
  cnb_plot <- ggplot(cnb_df, aes_string(x='Age', y=test, color='last_diagnosis')) +
    theme_linedraw() + geom_line(alpha=.5) +
    facet_grid(first_diagnosis ~ last_diagnosis) +
    scale_color_manual(values=c('deepskyblue2', 'darkorchid1', 'firebrick4')) +
    theme(legend.position = 'none') + labs(title=test)
  assign(paste0(test, '_plot'), cnb_plot)
}

# TO DO:
# https://stackoverflow.com/questions/31075407/plot-mixed-effects-model-in-ggplot
# Subtitle with Ns specific to the test


pdf(file='~/Documents/pncLongitudinalPsychosis/plots/longCog3x3.pdf', width=9, height=7)
ADT_plot
CPF_plot
CPT_plot
CPW_plot
ER40_plot
MEDF_plot
MPRAXIS_plot
NBACK_plot
PCET_plot
PLOT_plot
PMAT_plot
PVRT_plot
TAP_plot
VOLT_plot
dev.off()
