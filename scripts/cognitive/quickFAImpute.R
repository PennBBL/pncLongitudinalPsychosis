### This script does an exploratory factor analysis of the longitudinal CNB
### data, without accounting for repeated measures. As such, it is important
### that this ultimately be redone.
###
### Ellyn Butler
### August 11, 2020 - August 12, 2020 (Patch October 19, 2020)
### Fixed labels on March 2, 2021

set.seed(20)

library('dplyr') # Version 1.0.2
library('reshape2') # Version 1.4.4
library('ggplot2') # Version 3.3.2
library('ggpubr') # Version 0.4.0
library('psych') # Version 2.0.9
library('missForest') # Version 1.4
library('GPArotation') # Version 2014.11-1

clin_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/clinical/pnc_longitudinal_diagnosis_n749_20210112.csv')

# Recalculate ntimepoints (a lot of erroneous zeros)
for (bblid in unique(clin_df$bblid)) {
  clin_df[clin_df$bblid == bblid, 'ntimepoints'] <- length(clin_df[clin_df$bblid == bblid, 'ntimepoints'])
}

# Create first/last diagnoses df
clin_df <- clin_df[clin_df$timepoint == 't1' | clin_df$timepoint == paste0('t', clin_df$ntimepoints), ]
clin_df$timepoint <- recode(clin_df$timepoint, 't1'='t1', 't2'='tfinal2',
  't3'='tfinal2', 't4'='tfinal2', 't5'='tfinal2', 't6'='tfinal2')

clin_df$diagnosis <- recode(clin_df$diagnosis, 'psy'='PS')

clin_df <- reshape2::dcast(clin_df, bblid ~ timepoint, value.var='diagnosis')
clin_df$t1_tfinal <- paste(clin_df$t1, clin_df$tfinal2, sep='_')

clin_df$Diagnoses <- recode(clin_df$t1_tfinal, 'TD_TD'='TD-TD', 'TD_other'='TD-OP',
  'TD_PS'='TD-PS', 'other_TD'='OP-TD', 'other_other'='OP-OP', 'other_PS'='OP-PS',
  'PS_TD'='PS-TD', 'PS_other'='PS-OP', 'PS_PS'='PS-PS')

names(clin_df)[names(clin_df) == 't1'] <- 'first_diagnosis'
names(clin_df)[names(clin_df) == 'tfinal2'] <- 'last_diagnosis'
clin_df$first_diagnosis <- as.character(clin_df$first_diagnosis)
clin_df$last_diagnosis <- as.character(clin_df$last_diagnosis)
clin_df$first_diagnosis <- recode(clin_df$first_diagnosis, 'TD'='TD', 'other'='OP', 'PS'='PS')
clin_df$last_diagnosis <- recode(clin_df$last_diagnosis, 'TD'='TD', 'other'='OP', 'PS'='PS')
clin_df$t1_tfinal <- recode(clin_df$t1_tfinal, 'TD_TD'='TD_TD', 'TD_other'='TD_OP',
  'TD_PS'='TD_PS', 'other_TD'='OP_TD', 'other_other'='OP_OP', 'other_PS'='OP_PS',
  'PS_TD'='PS_TD', 'PS_other'='PS_OP', 'PS_PS'='PS_PS')

demo_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/demographics/baseline/n9498_demo_sex_race_ethnicity_dob.csv')
demo_df <- demo_df[, c('bblid', 'sex', 'race', 'ethnicity')]

cnb_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/cognitive/CNB_Longitudinal_Core_11February2020.csv')
cnb_df <- cnb_df[, c('bblid', 'Age', 'Test', 'ACC_raw', 'RT_raw')]
cnb_df <- cnb_df[cnb_df$bblid %in% clin_df$bblid,]

#### Calculate timepoint (ADDED OCTOBER 19, 2020) ####
bblids <- unique(cnb_df$bblid)

getTimepoint <- function(i) {
  bblid <- cnb_df[i, 'bblid']
  ages <- cnb_df[cnb_df$bblid == bblid, 'Age']
  ages <- sort(unique(ages))
  # Return timepoint
  which(ages == cnb_df[i, 'Age'])
}

cnb_df$timepoint <- sapply(1:nrow(cnb_df), getTimepoint)

#### Get rid of assessments done in the same month for the same person
cnb_df$bblid_Age_Test <- paste(cnb_df$bblid, cnb_df$Age, cnb_df$Test, sep='_')
cnb_df <- cnb_df[!duplicated(cnb_df$bblid_Age_Test),]








# Reverse code the RT data (want faster to be higher)
cnb_df$RT_raw <- -cnb_df$RT_raw

cnb_df1 <- dcast(cnb_df, bblid + Age + timepoint ~ Test, value.var='ACC_raw')
names(cnb_df1) <- c('bblid', 'Age', 'Timepoint',
  paste0(names(cnb_df1)[4:length(names(cnb_df1))], '_ACC'))
cnb_df2 <- dcast(cnb_df, bblid + Age + timepoint ~ Test, value.var='RT_raw')
names(cnb_df2) <- c('bblid', 'Age', 'Timepoint',
  paste0(names(cnb_df2)[4:length(names(cnb_df2))], '_RT'))
cnb_df <- merge(cnb_df1, cnb_df2)

cnb_df <- merge(cnb_df, demo_df)
cnb_df$sex <- recode(cnb_df$sex, `2`='Female', `1`='Male')

# Remove TAP_RT, and rename TAP_ACC to TAP_RT
cnb_df <- cnb_df[,!(names(cnb_df) %in% c('TAP_RT', 'MPRAXIS_ACC'))]
names(cnb_df)[names(cnb_df) == 'TAP_ACC'] <- 'TAP_RT'

tests_acc <- c('ADT', 'CPF', 'CPT', 'CPW', 'ER40', 'MEDF', 'NBACK', 'PCET',
  'PLOT', 'PMAT', 'PVRT', 'VOLT')
tests_rt <- c('ADT', 'CPF', 'CPT', 'CPW', 'ER40', 'MEDF', 'NBACK', 'PCET',
  'PLOT', 'PMAT', 'PVRT', 'VOLT', 'TAP', 'MPRAXIS')

cnb_df[, c(paste0(tests_acc, '_ACC'), paste0(tests_rt, '_RT'))] <- sapply(cnb_df[,
  c(paste0(tests_acc, '_ACC'), paste0(tests_rt, '_RT'))], scale)

getDiagnoses <- function(i) {
  bblid <- cnb_df[i, 'bblid']
  c(clin_df[clin_df$bblid == bblid, 'first_diagnosis'], clin_df[clin_df$bblid == bblid, 'last_diagnosis'], as.character(clin_df[clin_df$bblid == bblid, 't1_tfinal']))
}

cnb_df[,c('first_diagnosis', 'last_diagnosis', 't1_tfinal')] <- t(sapply(1:nrow(cnb_df), getDiagnoses))

cnb_df$first_diagnosis <- ordered(cnb_df$first_diagnosis, c('TD', 'OP', 'PS'))
cnb_df$last_diagnosis <- ordered(cnb_df$last_diagnosis, c('TD', 'OP', 'PS'))
cnb_df$t1_tfinal <- relevel(factor(cnb_df$t1_tfinal), ref='TD_TD')

tmp_df <- cnb_df

# Impute data
tmp_df[, c(paste0(tests_acc, '_ACC'), paste0(tests_rt, '_RT'))] <- missForest(tmp_df[,
  c(paste0(tests_acc, '_ACC'), paste0(tests_rt, '_RT'))])$ximp

# Calculate efficiency
calcEfficiency <- function(test) {
  rowMeans(tmp_df[, c(paste0(test, '_ACC'), paste0(test, '_RT'))])
}
tmp_df[, paste0(tests_acc, '_EFF')] <- sapply(tests_acc, calcEfficiency)
tmp_df[, paste0(tests_acc, '_EFF')] <- sapply(tmp_df[, paste0(tests_acc, '_EFF')], scale)

# Get factor scores and save loadings
types <- c('ACC', 'RT', 'EFF')
numfactors <- 1:4
tests_eff <- tests_acc

for (type in types) {
  tests <- get(paste0('tests_', tolower(type)))
  cnb_eigenvalues <- eigen(cor(tmp_df[, paste0(tests, '_', type)]))$values
  eigen_df <- data.frame(matrix(NA, nrow=length(cnb_eigenvalues), ncol=2))
  names(eigen_df) <- c("compnum", "eigen")
  eigen_df$compnum <- 1:length(tests)
  eigen_df$eigen <- cnb_eigenvalues

  eigen_plot <- ggplot(eigen_df, aes(x=compnum, y=eigen)) +
    geom_line(stat="identity") + geom_point() +  theme_minimal() +
    xlab("Component Number") + ylab("Eigenvalues of Components") +
    scale_y_continuous(limits=c(0, 6)) + ggtitle(paste0(type,': Scree Plot')) +
    theme(plot.title = element_text(size=12), axis.title = element_text(size=10),
      axis.text = element_text(size=6))
  assign(paste0(type, '_eigen_plot'), eigen_plot)

  for (numfs in numfactors) {
    fanal <- fa(tmp_df[, paste0(tests, '_', type)], nfactors=numfs, rotate='oblimin')
    if (numfs == 1) {
      loadings_df <- unclass(fanal$loadings)
      colnames(loadings_df) <- paste0('Soln', numfs, '_', colnames(loadings_df))
    } else {
      loadings_tmp_df <- unclass(fanal$loadings)
      colnames(loadings_tmp_df) <- paste0('Soln', numfs, '_', colnames(loadings_tmp_df))
      loadings_df <- cbind(loadings_df, loadings_tmp_df)
    }
    #factor_info <- factor.scores(tmp_df[, paste0(tests, '_', type)], fanal)
    tmp_df[, paste0(type, '_Soln', numfs, '_MR', 1:numfs)] <- fanal$scores
  }
  loadings_df <- loadings_df[, sort(colnames(loadings_df))]
  loadings_df <- as.data.frame(loadings_df)
  # Sort rows
  if (nrow(loadings_df) == 12) {
    loadings_df$desiredOrder <- c(10, 7, 4, 8, 11, 12, 5, 6, 1, 2, 3, 9)
    loadings_df <- loadings_df[order(loadings_df$desiredOrder), ]
    loadings_df$Domain <- c(rep('ComCog', 3), rep('Exec', 3), rep('Mem', 3), rep('SocCog', 3))
  } else {
    loadings_df$desiredOrder <- c(10, 7, 4, 8, 11, 12, 5, 6, 1, 2, 3, 9, 13, 14)
    loadings_df <- loadings_df[order(loadings_df$desiredOrder), ]
    loadings_df$Domain <- c(rep('ComCog', 3), rep('Exec', 3), rep('Mem', 3), rep('SocCog', 3), rep('SMSpeed', 2))
  }
  loadings_df <- loadings_df[, names(loadings_df) != 'desiredOrder']
  loadings_df$Test <- row.names(loadings_df)
  loadings_df <- loadings_df[, c(11, 12, 1:10)]
  assign(paste0(type, '_loadings_df'), loadings_df)
}

pdf(file='~/Documents/pncLongitudinalPsychosis/plots/eigenAccRtEff_impute.pdf', width=12, height=4)
ggarrange(ACC_eigen_plot, RT_eigen_plot, EFF_eigen_plot, ncol=3)
dev.off()

loadings_df <- rbind(ACC_loadings_df, RT_loadings_df, EFF_loadings_df)
loadings_df[, grep('MR', names(loadings_df), value=TRUE)] <- round(loadings_df[, grep('MR', names(loadings_df), value=TRUE)], digits=4)

write.csv(loadings_df, paste0('~/Documents/pncLongitudinalPsychosis/results/factorLoadings_impute_', Sys.Date(),'.csv'), row.names=FALSE)

# Select for only the columns you will be using in analysis
names(tmp_df)[names(tmp_df) == 'EFF_Soln4_MR1'] <- 'SocCog_EFF'
names(tmp_df)[names(tmp_df) == 'EFF_Soln4_MR2'] <- 'Exec_EFF'
names(tmp_df)[names(tmp_df) == 'EFF_Soln4_MR3'] <- 'Mem_EFF'
names(tmp_df)[names(tmp_df) == 'EFF_Soln4_MR4'] <- 'CompCog_EFF'

tmp_df <- tmp_df[, c('bblid', 'Age', 'race', 'sex', 'Timepoint', 'first_diagnosis', 'last_diagnosis',
  't1_tfinal', grep('_ACC', names(tmp_df), value=TRUE), grep('_RT', names(tmp_df), value=TRUE),
  grep('_EFF', names(tmp_df), value=TRUE))]

write.csv(tmp_df, paste0('~/Documents/pncLongitudinalPsychosis/data/cognitive/cnb_quickFA_impute_', Sys.Date(), '.csv'), row.names=FALSE)











#
