### This script does an exploratory factor analysis of the longitudinal CNB
### data, without accounting for repeated measures. As such, it is important
### that this ultimately be redone.
###
### Ellyn Butler
### August 11, 2020 - August 12, 2020

set.seed(20)

library('dplyr')
library('reshape2')
library('ggplot2')
library('ggpubr')
library('psych')

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

demo_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/demographics/baseline/n1601_demographics_go1_20161212.csv')

cnb_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/cognitive/CNB_Longitudinal_Core_11February2020.csv')
cnb_df <- cnb_df[, c('bblid', 'Age', 'timepoint', 'Test', 'ACC_raw', 'RT_raw')]
cnb_df <- cnb_df[cnb_df$bblid %in% clin_df$bblid,]

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
#cnb_df$t1_tfinal <- ordered(cnb_df$t1_tfinal, c('TD_TD', 'TD_OP', 'TD_PS',
#  'OP_TD', 'OP_OP', 'OP_PS', 'PS_TD', 'PS_OP', 'PS_PS'))
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
  assign(paste0(type, '_loadings_df'), loadings_df)
}

pdf(file='~/Documents/pncLongitudinalPsychosis/plots/eigenAccRtEff_impute_sex.pdf', width=12, height=4)
ggarrange(ACC_eigen_plot, RT_eigen_plot, EFF_eigen_plot, ncol=3)
dev.off()


loadings_df <- rbind(ACC_loadings_df, RT_loadings_df, EFF_loadings_df)
loadings_df <- round(loadings_df, digits=4)

write.csv(loadings_df, '~/Documents/pncLongitudinalPsychosis/results/factorLoadings_impute_sex.csv')

write.csv(tmp_df, '~/Documents/pncLongitudinalPsychosis/data/cnb_quickFA_impute_sex.csv', row.names=FALSE)
