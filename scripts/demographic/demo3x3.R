### This script makes plot of baseline demographic and trauma features by
### 3x3 trajectories
###
### Ellyn Butler
### July 23, 2020 - July 27, 2020

library('dplyr')
library('reshape2')
library('ggplot2')
library('ggpubr')
library('fastDummies')

clin_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/clinical/pnc_longitudinal_diagnosis_n752_202007.csv')
demo_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/demographics/baseline/n1601_demographics_go1_20161212.csv')
trauma_df <- read.csv('~/Documents/traumaInformant/data/PNC_GO1_GOASSESSDataArchiveNontext_DATA_2015-07-14_1157.csv')
trauma_df <- trauma_df[trauma_df$interview_type %in% c('AP', 'MP', 'YPI'),]
trauma_df <- trauma_df[,c('proband_bblid', grep('ptd', names(trauma_df), value=TRUE))]
names(trauma_df)[names(trauma_df) == 'proband_bblid'] <- 'bblid'

final_df <- merge(clin_df, demo_df, by='bblid')
final_df <- merge(final_df, trauma_df, by='bblid')
names(final_df)[names(final_df) == 't1'] <- 'first_diagnosis'
names(final_df)[names(final_df) == 'tfinal2'] <- 'last_diagnosis'
final_df$Female <- recode(final_df$sex, `2`=1, `1`=0)
final_df$White <- recode(final_df$race, `1`=1, .default=0)
final_df$first_diagnosis <- recode(final_df$first_diagnosis, 'other'='OP')
final_df$last_diagnosis <- recode(final_df$last_diagnosis, 'other'='OP')

ptdvars <- c(paste0('ptd00', 1:4), paste0('ptd00', 6:9))
final_df[, ptdvars] <- sapply(final_df[, ptdvars], na_if, y=9)

getPercent <- function(i, dataf) {
  first_diag <- as.character(dataf[i, 'first_diagnosis'])
  last_diag <- as.character(dataf[i, 'last_diagnosis'])
  feat <- as.character(dataf[i, 'Feature'])
  (nrow(final_df[final_df$first_diagnosis == first_diag & final_df$last_diagnosis ==
    last_diag & final_df[,feat] == 1  & !is.na(final_df[,feat]),])/nrow(final_df[final_df$first_diagnosis ==
    first_diag & final_df$last_diagnosis == last_diag & !is.na(final_df[,feat]),]))*100
}

subtit <- paste0('Ns: TD-TD=', nrow(final_df[final_df$t1_tfinal == 'TD_TD',]),
  ', TD-OP=', nrow(final_df[final_df$t1_tfinal == 'TD_other',]),
  ', TD-PS=', nrow(final_df[final_df$t1_tfinal == 'TD_PS',]),
  ', OP-TD=', nrow(final_df[final_df$t1_tfinal == 'other_TD',]),
  ', OP-OP=', nrow(final_df[final_df$t1_tfinal == 'other_other',]),
  ', OP-PS=', nrow(final_df[final_df$t1_tfinal == 'other_PS',]),
  ', PS-TD=', nrow(final_df[final_df$t1_tfinal == 'PS_TD',]),
  ', PS-OP=', nrow(final_df[final_df$t1_tfinal == 'PS_other',]),
  ', PS-PS=', nrow(final_df[final_df$t1_tfinal == 'PS_PS',]))

################################ Types ################################

final_type_df <- expand.grid(c('TD', 'OP', 'PS'), c('TD', 'OP', 'PS'),
  c('Female', 'White', paste0('ptd00', 1:4), paste0('ptd00', 6:9)))
names(final_type_df) <- c('first_diagnosis', 'last_diagnosis', 'Feature')

final_type_df$Percent <- sapply(1:nrow(final_type_df), getPercent, dataf=final_type_df)

final_type_df$Feature <- recode(final_type_df$Feature, 'ptd001'='disaster',
  'ptd002'='threat_close', 'ptd003'='physical', 'ptd004'='sexual',
  'ptd005'='rape', 'ptd006'='threat_weapon', 'ptd007'='accident',
  'ptd008'='witness', 'ptd009'='body') #something screwy is going on with rape variable
final_type_df$type <- recode(final_type_df$Feature, 'White'='demo', 'Female'='demo',
  'disaster'='ptd', 'threat_close'='ptd', 'physical'='ptd', 'sexual'='ptd',
  'rape'='ptd', 'threat_weapon'='ptd', 'accident'='ptd', 'witness'='ptd',
  'body'='ptd')

final_type_df$first_diagnosis <- recode(final_type_df$first_diagnosis,
  'PS'='PS - First Diagnosis', 'OP'='OP - First Diagnosis', 'TD'='TD - First Diagnosis')
final_type_df$first_diagnosis <- ordered(final_type_df$first_diagnosis,
  c('PS - First Diagnosis', 'OP - First Diagnosis', 'TD - First Diagnosis'))
final_type_df$last_diagnosis <- recode(final_type_df$last_diagnosis,
  'TD'='TD - Last Diagnosis', 'OP'='OP - Last Diagnosis', 'PS'='PS - Last Diagnosis')
final_type_df$last_diagnosis <- ordered(final_type_df$last_diagnosis,
  c('TD - Last Diagnosis', 'OP - Last Diagnosis', 'PS - Last Diagnosis'))

type_plot <- ggplot(final_type_df[final_type_df$Feature != 'rape',],
    aes(x=Feature, y=Percent, fill=type)) +
  theme_linedraw() + geom_bar(stat = 'identity') +
  facet_grid(first_diagnosis ~ last_diagnosis) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position='none') +
  coord_cartesian(ylim=c(0, 100)) +
  labs(title='Demographic and Traumatic Features by Diagnostic Bin', subtitle=subtit)


pdf(file='~/Documents/pncLongitudinalPsychosis/plots/demoTraumaType3x3.pdf', width=9, height=7)
type_plot
dev.off()




################################ Sum ################################

final_df$NumTypesTraumas <- rowSums(final_df[, ptdvars])
final_df$NumTypesTraumas <- recode(final_df$NumTypesTraumas, `0`=0, `1`=1, `2`=2,
  `3`=3, `4`=3, `5`=3, `6`=3, `7`=3)
final_df$NumTypesTraumas <- factor(final_df$NumTypesTraumas)
final_df <- fastDummies::dummy_cols(final_df, select_columns = 'NumTypesTraumas')

final_sum_df <- expand.grid(c('TD', 'OP', 'PS'), c('TD', 'OP', 'PS'),
  c('Female', 'White', paste0('NumTypesTraumas_', 0:3)))
names(final_sum_df) <- c('first_diagnosis', 'last_diagnosis', 'Feature')

final_sum_df$Percent <- sapply(1:nrow(final_sum_df), getPercent, dataf=final_sum_df)

final_sum_df$first_diagnosis <- recode(final_sum_df$first_diagnosis,
  'PS'='PS - First Diagnosis', 'OP'='OP - First Diagnosis', 'TD'='TD - First Diagnosis')
final_sum_df$first_diagnosis <- ordered(final_sum_df$first_diagnosis,
  c('PS - First Diagnosis', 'OP - First Diagnosis', 'TD - First Diagnosis'))
final_sum_df$last_diagnosis <- recode(final_sum_df$last_diagnosis,
  'TD'='TD - Last Diagnosis', 'OP'='OP - Last Diagnosis', 'PS'='PS - Last Diagnosis')
final_sum_df$last_diagnosis <- ordered(final_sum_df$last_diagnosis,
  c('TD - Last Diagnosis', 'OP - Last Diagnosis', 'PS - Last Diagnosis'))

final_sum_df$type <- recode(final_sum_df$Feature, 'White'='demo', 'Female'='demo',
    'NumTypesTraumas_0'='ptd', 'NumTypesTraumas_1'='ptd', 'NumTypesTraumas_2'='ptd',
    'NumTypesTraumas_3'='ptd')

sum_plot <- ggplot(final_sum_df, aes(x=Feature, y=Percent, fill=type)) +
  theme_linedraw() + geom_bar(stat = 'identity') +
  facet_grid(first_diagnosis ~ last_diagnosis) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position='none') +
  coord_cartesian(ylim=c(0, 100)) +
  labs(title='Demographic and Traumatic Features by Diagnostic Bin', subtitle=subtit)


pdf(file='~/Documents/pncLongitudinalPsychosis/plots/demoTraumaSum3x3.pdf', width=9, height=7)
sum_plot
dev.off()

################################ Assault ################################

# Assault, Threat, Other, None
final_df$assault <- rowSums(final_df[, c('ptd003', 'ptd004')])
final_df$assault <- recode(final_df$assault, `0`=0, `1`=1, `2`=1)
final_df$threat <- rowSums(final_df[, c('ptd002', 'ptd006')])
final_df$threat <- recode(final_df$threat, `0`=0, `1`=1, `2`=1)
final_df$other <- rowSums(final_df[, c('ptd001', 'ptd007', 'ptd008', 'ptd009')])
final_df$other <- recode(final_df$other, `0`=0, `1`=1, `2`=1, `3`=1, `4`=1)

final_assault_df <- expand.grid(c('TD', 'OP', 'PS'), c('TD', 'OP', 'PS'),
  c('Female', 'White', 'assault', 'threat', 'other'))
names(final_assault_df) <- c('first_diagnosis', 'last_diagnosis', 'Feature')

final_assault_df$Percent <- sapply(1:nrow(final_assault_df), getPercent, dataf=final_assault_df)

final_assault_df$first_diagnosis <- recode(final_assault_df$first_diagnosis,
  'PS'='PS - First Diagnosis', 'OP'='OP - First Diagnosis', 'TD'='TD - First Diagnosis')
final_assault_df$first_diagnosis <- ordered(final_assault_df$first_diagnosis,
  c('PS - First Diagnosis', 'OP - First Diagnosis', 'TD - First Diagnosis'))
final_assault_df$last_diagnosis <- recode(final_assault_df$last_diagnosis,
  'TD'='TD - Last Diagnosis', 'OP'='OP - Last Diagnosis', 'PS'='PS - Last Diagnosis')
final_assault_df$last_diagnosis <- ordered(final_assault_df$last_diagnosis,
  c('TD - Last Diagnosis', 'OP - Last Diagnosis', 'PS - Last Diagnosis'))

final_assault_df$type <- recode(final_assault_df$Feature, 'White'='demo', 'Female'='demo',
    'assault'='ptd', 'threat'='ptd', 'other'='ptd')

assault_plot <- ggplot(final_assault_df, aes(x=Feature, y=Percent, fill=type)) +
  theme_linedraw() + geom_bar(stat = 'identity') +
  facet_grid(first_diagnosis ~ last_diagnosis) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position='none') +
  coord_cartesian(ylim=c(0, 100)) +
  labs(title='Demographic and Traumatic Features by Diagnostic Bin', subtitle=subtit)


pdf(file='~/Documents/pncLongitudinalPsychosis/plots/demoTraumaAssault3x3.pdf', width=9, height=7)
assault_plot
dev.off()
