### This script plots fits and outputs a model summary table for the Freesurfer
### imaging values.
### NOTE TO FUTURE DATA ANALYSTS: This will need to be redone with ComBat + GAMM.
### Reach out to Dr. Joanne Beer for help.
###
### Ellyn Butler
### June 1, 2021 - June 2, 2021

set.seed(20)

library('dplyr') # V 1.0.2
library('reshape2') # V 1.4.4
library('ggplot2') # V 3.3.2
library('ggpubr') # V 0.4.0
library('lme4') # V 1.1-26
library('gamm4') # V 0.2-6
library('stringr') # V 1.4.0
library('MASS') # V 7.3-53
library('broom') # V 0.7.2
library('sjPlot') # V 2.8.6
library('cowplot') # V 1.1.1
library('pbkrtest') # V 0.5-0.1

# Declare functions
getAssessmentNumber <- function(i, dataf) {
  bblid <- dataf[i, 'bblid']
  ages <- dataf[dataf$bblid == bblid, 'age']
  ages <- sort(unique(ages))
  which(ages == dataf[i, 'age'])
}

getUpperLowerCI <- function(i) {
  sorted_vec <- unname(sort(fits[i,]))
  lower <- sorted_vec[round(.025*length(sorted_vec))]
  upper <- sorted_vec[round(.975*length(sorted_vec))]
  c(lower, upper)
}

# Read in data
area_lh_df <- read.csv('~/Documents/ExtraLong/data/freesurferCrossSectional/tabulated/lh_DKTatlas_area_2020-11-09.csv')
cort_lh_df <- read.csv('~/Documents/ExtraLong/data/freesurferCrossSectional/tabulated/lh_DKTatlas_thickness_2020-11-09.csv')
area_rh_df <- read.csv('~/Documents/ExtraLong/data/freesurferCrossSectional/tabulated/rh_DKTatlas_area_2020-11-09.csv')
cort_rh_df <- read.csv('~/Documents/ExtraLong/data/freesurferCrossSectional/tabulated/rh_DKTatlas_thickness_2020-11-09.csv')
area_df <- merge(area_lh_df, area_rh_df)
cort_df <- merge(cort_lh_df, cort_rh_df)
img_df <- merge(area_df, cort_df)
qual_df <- read.csv('~/Documents/ExtraLong/data/qualityAssessment/antssstExclude.csv')
img_df <- merge(qual_df, img_df)
img_df <- img_df[img_df$antssstExclude == FALSE, ]
row.names(img_df) <- 1:nrow(img_df)


diag_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/clinical/pnc_longitudinal_diagnosis_n749_20210112.csv')
demo_df <- read.csv('~/Documents/ExtraLong/data/demographicsClinical/scanid_to_seslabel_demo_20200531.csv')

# Recalculate ntimepoints (a lot of erroneous zeros)
for (bblid in unique(diag_df$bblid)) {
  diag_df[diag_df$bblid == bblid, 'ntimepoints'] <- length(diag_df[diag_df$bblid == bblid, 'ntimepoints'])
}

# Create first/last diagnoses df
diag_df <- diag_df[diag_df$timepoint == 't1' | diag_df$timepoint == paste0('t', diag_df$ntimepoints), ]
diag_df$timepoint <- recode(diag_df$timepoint, 't1'='t1', 't2'='tfinal2',
  't3'='tfinal2', 't4'='tfinal2', 't5'='tfinal2', 't6'='tfinal2')

diag_df$diagnosis <- recode(diag_df$diagnosis, 'psy'='PS')

diag_df <- reshape2::dcast(diag_df, bblid ~ timepoint, value.var='diagnosis')
diag_df$t1_tfinal <- paste(diag_df$t1, diag_df$tfinal2, sep='_')

diag_df$Diagnoses <- recode(diag_df$t1_tfinal, 'TD_TD'='TD-TD', 'TD_other'='TD-OP',
  'TD_PS'='TD-PS', 'other_TD'='OP-TD', 'other_other'='OP-OP', 'other_PS'='OP-PS',
  'PS_TD'='PS-TD', 'PS_other'='PS-OP', 'PS_PS'='PS-PS')

# Name/relevel variables as necessary
diag_df$Diagnoses <- factor(diag_df$Diagnoses)
diag_df$Diagnoses <- relevel(diag_df$Diagnoses, 'TD-TD')

diag_df$oDiagnoses <- ordered(diag_df$Diagnoses, c('TD-TD', 'OP-OP', 'OP-PS',
  'OP-TD', 'PS-OP', 'PS-PS', 'PS-TD', 'TD-OP', 'TD-PS'))

img_df <- merge(img_df, diag_df, by='bblid')
img_df <- merge(img_df, demo_df)

img_df$White <- recode(img_df$race, `1`=TRUE, `2`=FALSE, `3`=FALSE, `4`=FALSE, `5`=FALSE)
img_df$Male <- recode(img_df$sex, `1`=TRUE, `2`=FALSE)

names(img_df)[names(img_df) == 'scanage_years'] <- 'Age'



################################# Plot regions #################################

plotcols <- c(grep('_thickness', names(img_df), value=TRUE), grep('_area', names(img_df), value=TRUE))

region <- sapply(strsplit(plotcols, '_'), `[[`, 2)
region <- region[!(region %in% c('MeanThickness', 'WhiteSurfArea'))]
regionlobe_df <- data.frame(region=unique(region), lobe=c('cingulate', 'frontal',
  'occipital', 'temporal', 'temporal', 'parietal', 'temporal', 'cingulate',
  'occipital', 'frontal', 'occipital', 'frontal', 'temporal', 'temporal',
  'frontal', 'frontal', 'frontal', 'frontal', 'occipital', 'parietal',
  'cingulate', 'frontal', 'parietal', 'cingulate', 'frontal', 'frontal',
  'parietal', 'temporal', 'parietal', 'temporal', 'other'))

diags <- as.character(unique(img_df$Diagnoses)[!(unique(img_df$Diagnoses) %in% c('TD-TD', 'PS-PS'))])

i=1
for (Value in plotcols) {
  img_Value_df <- img_df
  img_Value_df <- img_Value_df[!is.na(img_Value_df[, Value]), ]
  row.names(img_Value_df) <- 1:nrow(img_Value_df)

  mod1b <- gamm4(as.formula(paste(Value, "~ Diagnoses + s(Age, k=4, bs='cr') +
    s(Age, by=oDiagnoses, k=4, bs='cr')")), data=img_Value_df, random=~(1|bblid), REML=TRUE)

  mod2b <- gamm4(as.formula(paste(Value, "~ Male + White + Diagnoses + s(Age, k=4, bs='cr') +
    s(Age, by=oDiagnoses, k=4, bs='cr')")), data=img_Value_df, random=~(1|bblid), REML=TRUE)

  print(tab_model(mod1b$gam, file=paste0('~/Documents/pncLongitudinalPsychosis/results/imaging/table_', Value, '_freesurfer.html')))
  assign(paste0(Value, '_model'), mod1b$gam)
  print(tab_model(mod2b$gam, file=paste0('~/Documents/pncLongitudinalPsychosis/results/imaging/tableSensitivity_', Value, '_freesurfer.html')))
  assign(paste0(Value, 'Sensitivity_model'), mod2b$gam)

  lp <- predict(mod1b$gam, newdata=img_Value_df, type='lpmatrix')
  coefs <- coef(mod1b$gam)
  vc <- vcov(mod1b$gam)

  sim <- mvrnorm(1000, mu = coefs, Sigma = vc)
  fits <- lp %*% t(sim)

  cis <- t(sapply(1:nrow(fits), getUpperLowerCI))
  cis <- data.frame(cis)
  names(cis) <- c('LCI', 'UCI')
  img_Value_df <- cbind(img_Value_df, cis)

  img_Value_df$predgamm <- predict(mod1b$gam)

  for (group in diags) {
    subtit <- paste0('Sessions: TD-TD=', nrow(img_Value_df[img_Value_df$Diagnoses == 'TD-TD',]),
      ', ', group, '=', nrow(img_Value_df[img_Value_df$Diagnoses == group,]),
      ', PS-PS=', nrow(img_Value_df[img_Value_df$Diagnoses == 'PS-PS',]))
    group_df <- img_Value_df[img_Value_df$Diagnoses %in% c('TD-TD', group, 'PS-PS'), ]
    group_df$Diagnoses <- ordered(group_df$Diagnoses, c('TD-TD', group, 'PS-PS'))

    intervals_df <- reshape2::melt(group_df, c('bblid', 'Age', 'Diagnoses'), c('LCI', 'UCI'))
    img_plot <- ggplot(group_df, aes_string(x='Age', y=Value, color='Diagnoses')) +
      theme_linedraw() + xlim(10, 30) +
      labs(title=Value, subtitle=subtit, y='Value (95% CI)', color='Diagnoses') +
      scale_color_manual(values=c('green3', 'goldenrod2', 'red')) +
      theme(legend.position = 'bottom', plot.title=element_text(size=12, face="bold"),
        plot.subtitle=element_text(size=7)) +
      geom_line(aes(y=predgamm), size=1) +
      geom_ribbon(data=img_Value_df[img_Value_df$Diagnoses == 'TD-TD',],
        aes(ymin=LCI, ymax=UCI), alpha=.2, fill='green3', show.legend=FALSE) +
      geom_ribbon(data=img_Value_df[img_Value_df$Diagnoses == group,],
        aes(ymin=LCI, ymax=UCI), alpha=.2, fill='goldenrod2', show.legend=FALSE) +
      geom_ribbon(data=img_Value_df[img_Value_df$Diagnoses == 'PS-PS',],
        aes(ymin=LCI, ymax=UCI), alpha=.2, fill='red', show.legend=FALSE)

    assign(paste0(Value, '_', gsub('-', '_', group), '_plot'), img_plot)

    pdf(file=paste0('~/Documents/pncLongitudinalPsychosis/plots/imaging/', Value, '_', group, '_freesurfer.pdf'), width=4, height=4)
    print(img_plot)
    dev.off()
  }
  i=i+1
}


# Create tables for each of the lobes
for (lobe in unique(regionlobe_df$lobe)) {
  for (modal in c('area', 'thickness')) {
    for (hemi in c('rh', 'lh')) {
      lobe_regions <- regionlobe_df[regionlobe_df$lobe == lobe, 'region']
      model_list <- list()
      for (i in 1:length(lobe_regions)) {
        reg <- lobe_regions[i]
        model_list[[i]] <- get(paste(hemi, reg, modal, 'model', sep='_'))
      }
      print(tab_model(model_list, p.adjust='fdr',
        file=paste0('~/Documents/pncLongitudinalPsychosis/results/imaging/', lobe, '_', modal, '_', hemi, '_freesurfer.html')))
    }
  }
}






#
