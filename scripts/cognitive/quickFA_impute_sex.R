### This script does an exploratory factor analysis of the longitudinal CNB
### data, without accounting for repeated measures. In addition, it imputes
### missing data. It plots the data by sex, tests for sex differences within
### diagnostic groups, and includes MPRAXIS (speed) and TAP (speed) in the
### factor analysis.
###
### Ellyn Butler
### August 24, 2020

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




###################### Plot factor scores ######################

plotcols <- c(paste0('ACC_Soln3_MR', 1:3), paste0('RT_Soln3_MR', 1:3),
  paste0('EFF_Soln4_MR', 1:4))

for (test in plotcols) {
  cnb_test_df <- tmp_df
  names(cnb_test_df)[names(cnb_test_df) == test] <- 'test'

  cnb_test_df$sex <- factor(cnb_test_df$sex)
  cnb_test_df$t1_tfinal_factor <- factor(cnb_test_df$t1_tfinal)

  #cnb_test_df$sex_t1_tfinal_factor <- paste0(cnb_test_df$sex, cnb_test_df$t1_tfinal_factor)

  #mod1b <- gamm4(test ~ s(Age, k=10, bs='cr') + s(Age, by=sex, k=10, bs='cr') +
  #  s(Age, by=t1_tfinal_factor, k=10, bs='cr') +
  #  s(Age, by=sex*t1_tfinal_factor, k=10, bs='cr'),
  #  data=cnb_test_df, random=~(1|bblid), REML=TRUE)
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
    scale_color_manual(values=c('blue', 'red')) +
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
