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

tests <- c('ADT', 'CPF', 'CPT', 'CPW', 'ER40', 'MEDF', 'NBACK', 'PCET',
  'PLOT', 'PMAT', 'PVRT', 'VOLT') # August 12, 2020: Got rid of TAP

cnb_df[, tests] <- sapply(cnb_df[, c(paste0(tests, '_ACC'), paste0(tests, '_RT'))], scale)
tmp_df <- cnb_df
tmp_df <- tmp_df[!is.na(tmp_df$ADT) & !is.na(tmp_df$CPF) & !is.na(tmp_df$CPT) &
  !is.na(tmp_df$CPW) & !is.na(tmp_df$ER40) & !is.na(tmp_df$MEDF) &
  !is.na(tmp_df$NBACK) & !is.na(tmp_df$PCET) & !is.na(tmp_df$PLOT) &
  !is.na(tmp_df$PMAT) & !is.na(tmp_df$PVRT) & !is.na(tmp_df$VOLT),]
  # August 11, 2020: EEK. Lose 789 if get rid of rows with any NA

#VSS.scree(tmp_df[,tests])

# Calculate efficiency
calcEfficiency <- function(test) {
  rowMeans(tmp_df[, c(paste0(test, '_ACC'), paste0(test, '_RT'))])
}
tmp_df[, paste0(tests, '_EFF')] <- sapply(tests, calcEfficiency)




# Get factor scores and save loadings
types <- c('ACC', 'RT', 'EFF')
numfactors <- 1:4

for (type in types) {
  cnb_eigenvalues <- eigen(cor(tmp_df[, paste0(tests, '_', type)]))$values
  eigen_df <- data.frame(matrix(NA, nrow=length(cnb_eigenvalues), ncol=2))
  names(eigen_df) <- c("compnum", "eigen")
  eigen_df$compnum <- 1:length(tests)
  eigen_df$eigen <- cnb_eigenvalues

  eigen_plot <- ggplot(eigen_df, aes(x=compnum, y=eigen)) +
    geom_line(stat="identity") + geom_point() +  theme_minimal() +
    xlab("Component Number") + ylab("Eigenvalues of Components") +
    scale_y_continuous(limits=c(0, 5)) + ggtitle(paste0(type,': Scree Plot')) +
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
    factor_info <- factor.scores(tmp_df[, paste0(tests, '_', type)], fanal)
    tmp_df[, paste0('MR', 1:numfs)] <- factor_info$scores
  }
  loadings_df <- loadings_df[, sort(colnames(loadings_df))]
  assign(paste0(type, '_loadings_df'), loadings_df)
}

pdf(file='~/Documents/pncLongitudinalPsychosis/plots/eigenAccRtEff.pdf', width=12, height=4)
ggarrange(ACC_eigen_plot, RT_eigen_plot, EFF_eigen_plot, ncol=3)
dev.off()


loadings_df <- rbind(ACC_loadings_df, RT_loadings_df, EFF_loadings_df)
loadings_df <- round(loadings_df, digits=4)

write.csv(loadings_df, '~/Documents/pncLongitudinalPsychosis/results/factorLoadings.csv')








write.csv(cnb_df, paste0('~/Documents/pncLongitudinalPsychosis/data/cognitive/factorscores_', Sys.Date(),'.csv'))
