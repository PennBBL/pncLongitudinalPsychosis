### This script makes plot of baseline demographic and trauma features by
### 3x3 trajectories
###
### Ellyn Butler
### July 23, 2020 - August 11, 2020

library('dplyr')
library('reshape2')
library('ggplot2')
library('ggpubr')
library('fastDummies')
library('broom')
library('stringr')

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
final_df$t1_tfinal <- recode(final_df$t1_tfinal, 'other_other'='OP_OP',
  'other_TD'='OP_TD', 'other_PS'='OP_PS', 'TD_other'='TD_OP', 'PS_other'='PS_OP')
final_df <- within(final_df, t1_tfinal <- relevel(t1_tfinal, ref='TD_TD'))

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
    ', TD-OP=', nrow(final_df[final_df$t1_tfinal == 'TD_OP',]),
    ', TD-PS=', nrow(final_df[final_df$t1_tfinal == 'TD_PS',]),
    ',\nOP-TD=', nrow(final_df[final_df$t1_tfinal == 'OP_TD',]),
    ', OP-OP=', nrow(final_df[final_df$t1_tfinal == 'OP_OP',]),
    ', OP-PS=', nrow(final_df[final_df$t1_tfinal == 'OP_PS',]),
    ',\nPS-TD=', nrow(final_df[final_df$t1_tfinal == 'PS_TD',]),
    ', PS-OP=', nrow(final_df[final_df$t1_tfinal == 'PS_OP',]),
    ', PS-PS=', nrow(final_df[final_df$t1_tfinal == 'PS_PS',]))

################################ Types ################################

final_type_df <- expand.grid(c('TD', 'OP', 'PS'), c('TD', 'OP', 'PS'),
  c('Female', 'White', ptdvars))
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
final_df$NumTypesTraumas2 <- final_df$NumTypesTraumas + .1
mod_sum <- glm(NumTypesTraumas2 ~ I(t1_tfinal), data=final_df, family=Gamma())
model_info <- tidy(mod_sum) %>% filter(str_detect(term, "t1_tfinal"))
model_info <- model_info[model_info$p.value < .05,]
diagcats <- gsub('factor', '', model_info$term)
diagcats <- gsub('I\\(t1_tfinal\\)', '', diagcats)
ann_text <- data.frame(t1_tfinal=diagcats, lab = '*', Feature='2+')
int <- strsplit(as.character(ann_text$t1_tfinal), '_')
ann_text$first_diagnosis <- sapply(int, `[[`, 1)
ann_text$first_diagnosis <- paste(ann_text$first_diagnosis, '- First Diagnosis')
ann_text$last_diagnosis <- sapply(int, `[[`, 2)
ann_text$last_diagnosis <- paste(ann_text$last_diagnosis, '- Last Diagnosis')
ann_text$first_diagnosis <- recode(ann_text$first_diagnosis,
  'TD'='TD - First Diagnosis', 'OP'='OP - First Diagnosis', 'PS'='PS - First Diagnosis')
ann_text$first_diagnosis <- ordered(ann_text$first_diagnosis,
  c('TD - First Diagnosis', 'OP - First Diagnosis', 'PS - First Diagnosis'))
ann_text$last_diagnosis <- recode(ann_text$last_diagnosis,
  'TD'='TD - Last Diagnosis', 'OP'='OP - Last Diagnosis', 'PS'='PS - Last Diagnosis')
ann_text$last_diagnosis <- ordered(ann_text$last_diagnosis,
  c('TD - Last Diagnosis', 'OP - Last Diagnosis', 'PS - Last Diagnosis'))


final_df$NumTypesTraumas <- recode(final_df$NumTypesTraumas, `0`=0, `1`=1, `2`=2,
  `3`=2, `4`=2, `5`=2, `6`=2, `7`=2)
final_df$NumTypesTraumas <- factor(final_df$NumTypesTraumas)
final_df <- fastDummies::dummy_cols(final_df, select_columns = 'NumTypesTraumas')

final_sum_df <- expand.grid(c('TD', 'OP', 'PS'), c('TD', 'OP', 'PS'),
  c(paste0('NumTypesTraumas_', 0:2)))
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

final_sum_df$Feature <- as.character(final_sum_df$Feature)
sum_plot <- ggplot(final_sum_df, aes(x=Feature, y=Percent, fill=Feature)) +
  theme_linedraw() + geom_bar(stat = 'identity') +
  facet_grid(first_diagnosis ~ last_diagnosis) +
  theme(legend.position='none') + coord_cartesian(ylim=c(0, 100)) +
  scale_x_discrete(breaks=c('NumTypesTraumas_0', 'NumTypesTraumas_1',
    'NumTypesTraumas_2'), labels=c('0', '1', '2+')) +
  scale_fill_manual(values=c('white', 'green3', 'red', 'red')) +
  labs(title='Number of Types of Traumas', subtitle=subtit) +
  geom_text(data=ann_text, y=80, label="*", size=15) # August 11, 2020: Can't get colors and asterisks to work at the same time


pdf(file='~/Documents/pncLongitudinalPsychosis/plots/traumaSum3x3.pdf', width=5, height=5.5)
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
final_df$none <- rowSums(final_df[, ptdvars])
final_df$none <- recode(final_df$other, `0`=1, `1`=0, `2`=0, `3`=0, `4`=0, `5`=0, `6`=0, `7`=0)

# Conduct proportion tests
ann_text <- expand.grid(c('TD', 'OP', 'PS'), c('TD', 'OP', 'PS'),
  c('none','assault', 'threat', 'other'))
names(ann_text) <- c('first_diagnosis', 'last_diagnosis', 'Feature')
ann_text$t1_tfinal <- paste(ann_text$first_diagnosis, ann_text$last_diagnosis, sep='_')
ann_text <- ann_text[ann_text$t1_tfinal != 'TD_TD',]
row.names(ann_text) <- 1:nrow(ann_text)

testPropTDTD <- function(i) {
  feat <- as.character(ann_text[i, 'Feature'])
  group <- as.character(ann_text[i, 't1_tfinal'])
  N_tdtd <- nrow(final_df[final_df$t1_tfinal == 'TD_TD' & !is.na(final_df[,feat]),])
  X_tdtd <- nrow(final_df[final_df$t1_tfinal == 'TD_TD' & final_df[,feat] == 1
    & !is.na(final_df[,feat]),])
  N_group <- nrow(final_df[final_df$t1_tfinal == group & !is.na(final_df[,feat]),])
  X_group <- nrow(final_df[final_df$t1_tfinal == group & final_df[,feat] == 1
    & !is.na(final_df[,feat]),])
  prop.test(c(X_tdtd, X_group), c(N_tdtd, N_group))$p.value
}
ann_text$p.value <- sapply(1:nrow(ann_text), testPropTDTD)
ann_text$p.value <- p.adjust(ann_text$p.value, method='hochberg') # for independent tests
ann_text <- ann_text[ann_text$p.value < .05,]
ann_text$first_diagnosis <- recode(ann_text$first_diagnosis,
  'TD'='TD - First Diagnosis', 'OP'='OP - First Diagnosis', 'PS'='PS - First Diagnosis')
ann_text$first_diagnosis <- ordered(ann_text$first_diagnosis,
  c('TD - First Diagnosis', 'OP - First Diagnosis', 'PS - First Diagnosis'))
ann_text$last_diagnosis <- recode(ann_text$last_diagnosis,
  'TD'='TD - Last Diagnosis', 'OP'='OP - Last Diagnosis', 'PS'='PS - Last Diagnosis')
ann_text$last_diagnosis <- ordered(ann_text$last_diagnosis,
  c('TD - Last Diagnosis', 'OP - Last Diagnosis', 'PS - Last Diagnosis'))
ann_text$lab <- '*'






final_assault_df <- expand.grid(c('TD', 'OP', 'PS'), c('TD', 'OP', 'PS'),
  c('none','assault', 'threat', 'other'))
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

assault_plot <- ggplot(final_assault_df, aes(x=Feature, y=Percent, fill=Feature)) +
  theme_linedraw() + geom_bar(stat = 'identity') +
  facet_grid(first_diagnosis ~ last_diagnosis) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position='none') +
  coord_cartesian(ylim=c(0, 100)) + scale_fill_manual(values=c('green3', 'red', 'red', 'red')) +
  labs(title='Trauma Types by Diagnostic Bin', subtitle=subtit) +
  geom_text(data=ann_text, y=80, label="*", size=15)


pdf(file='~/Documents/pncLongitudinalPsychosis/plots/trauma3x3.pdf', width=12, height=7)
ggarrange(sum_plot, assault_plot, ncol=2)
dev.off()
