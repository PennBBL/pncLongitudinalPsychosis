### This script plots age at worst experience by longitudinal clinical labels
###
### Ellyn Butler
### August 12, 2020 - August 24, 2020

library('dplyr')
library('reshape2')
library('ggplot2')
library('ggpubr')

clin_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/clinical/pnc_longitudinal_diagnosis_n752_202007.csv')

trauma_df <- read.csv('~/Documents/traumaInformant/data/PNC_GO1_GOASSESSDataArchiveNontext_DATA_2015-07-14_1157.csv')
trauma_df <- trauma_df[trauma_df$interview_type %in% c('AP', 'MP', 'YPI'),]
trauma_df <- trauma_df[,c('proband_bblid', grep('ptd', names(trauma_df), value=TRUE))]
names(trauma_df)[names(trauma_df) == 'proband_bblid'] <- 'bblid'

final_df <- merge(clin_df, trauma_df, by='bblid')
final_df <- final_df[(final_df$ptd001 == 1 | final_df$ptd002 == 1 |
  final_df$ptd003 == 1 | final_df$ptd004 == 1 | final_df$ptd006 == 1  |
  final_df$ptd007 == 1 | final_df$ptd008 == 1 | final_df$ptd009 == 1) &
  !is.na(final_df$ptd020),] ## NAs generated here

final_df$t1_tfinal <- recode(final_df$t1_tfinal, 'other_other'='OP_OP',
  'other_TD'='OP_TD', 'other_PS'='OP_PS', 'TD_other'='TD_OP', 'PS_other'='PS_OP')
final_df <- within(final_df, t1_tfinal <- relevel(t1_tfinal, ref='TD_TD'))
names(final_df)[names(final_df) == 't1'] <- 'first_diagnosis'
names(final_df)[names(final_df) == 'tfinal2'] <- 'last_diagnosis'
final_df$first_diagnosis <- recode(final_df$first_diagnosis, 'other'='OP')
final_df$last_diagnosis <- recode(final_df$last_diagnosis, 'other'='OP')

###############################################################################
# PTD020: "When did (insert worst event name) occur? Age"

final_df$first_diagnosis <- recode(final_df$first_diagnosis,
  'PS'='PS - First Diagnosis', 'OP'='OP - First Diagnosis', 'TD'='TD - First Diagnosis')
final_df$last_diagnosis <- recode(final_df$last_diagnosis,
  'TD'='TD - Last Diagnosis', 'OP'='OP - Last Diagnosis', 'PS'='PS - Last Diagnosis')

subtit <- paste0('Ns: TD-TD=', nrow(final_df[final_df$t1_tfinal == 'TD_TD',]),
    ', TD-OP=', nrow(final_df[final_df$t1_tfinal == 'TD_OP',]),
    ', TD-PS=', nrow(final_df[final_df$t1_tfinal == 'TD_PS',]),
    ',\nOP-TD=', nrow(final_df[final_df$t1_tfinal == 'OP_TD',]),
    ', OP-OP=', nrow(final_df[final_df$t1_tfinal == 'OP_OP',]),
    ', OP-PS=', nrow(final_df[final_df$t1_tfinal == 'OP_PS',]),
    ',\nPS-TD=', nrow(final_df[final_df$t1_tfinal == 'PS_TD',]),
    ', PS-OP=', nrow(final_df[final_df$t1_tfinal == 'PS_OP',]),
    ', PS-PS=', nrow(final_df[final_df$t1_tfinal == 'PS_PS',]))

vline_df <- expand.grid(unique(as.character(final_df$first_diagnosis)),
  unique(as.character(final_df$last_diagnosis)))
names(vline_df) <- c('first_diagnosis', 'last_diagnosis')
meansGroups <- function(i) {
  first_diag <- vline_df[i, 'first_diagnosis']
  last_diag <- vline_df[i, 'last_diagnosis']
  mean(final_df[final_df$first_diagnosis == first_diag & final_df$last_diagnosis == last_diag, 'ptd020'])
}
vline_df$Means <- sapply(1:nrow(vline_df), meansGroups)

final_df$first_diagnosis <- ordered(final_df$first_diagnosis,
  c('PS - First Diagnosis', 'OP - First Diagnosis', 'TD - First Diagnosis'))
final_df$last_diagnosis <- ordered(final_df$last_diagnosis,
  c('TD - Last Diagnosis', 'OP - Last Diagnosis', 'PS - Last Diagnosis'))

worst_plot <- ggplot(final_df, aes(x=ptd020)) + theme_linedraw() +
  geom_histogram(bins=20) + facet_grid(first_diagnosis ~ last_diagnosis) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position='none') +
  labs(title='Age at Worst Traumatic Experience', subtitle=subtit, x='Age at Worst') +
  geom_vline(data = vline_df, aes(xintercept=Means), color="blue", linetype="dashed", size=1)

pdf(file='~/Documents/pncLongitudinalPsychosis/plots/ageAtWorst.pdf', width=5, height=6)
worst_plot
dev.off()


# Plot means - per facet
# Figure out where NA diagnoses are coming from
# Order diagnoses
