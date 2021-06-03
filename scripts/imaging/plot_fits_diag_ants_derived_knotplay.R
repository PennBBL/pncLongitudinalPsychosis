### This script plots fits and outputs a model summary table for the ANTs imaging
### values. I played with the number and placement of knots to see if that changed
### the lobular volume curves. I did this because the curve does not match the
### lower age range, and Ruben wanted to make sure that what we were finding
### wasn't an artifact of the method.
### NOTE TO FUTURE DATA ANALYSTS: This will need to be redone with ComBat + GAMM.
### Reach out to Dr. Joanne Beer for help.
###
### Ellyn Butler
### June 2, 2021

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

############################### Declare functions ###############################

getUpperLowerCI <- function(i) {
  sorted_vec <- unname(sort(fits[i,]))
  lower <- sorted_vec[round(.025*length(sorted_vec))]
  upper <- sorted_vec[round(.975*length(sorted_vec))]
  c(lower, upper)
}

getMoreRegions <- function(lobe) {
  regions <- regionlobe_df[regionlobe_df$lobe == lobe, 'region']
  regions_vol_rh <- paste0('vol_rh_', regions)
  regions_ct_rh <- paste0('ct_rh_', regions)
  regions_gmd_rh <- paste0('gmd_rh_', regions)

  getAveDiffVol <- function(region_rh) {
    region_lh <- gsub('_rh_', '_lh_', region_rh)
    tmp <- data.frame(
      ave=rowMeans(img_df[, c(region_rh, region_lh)]),
      diff=img_df[, region_rh] - img_df[, region_lh]
    )
    names(tmp) <- c(gsub('_rh_', '_ave_', region_rh), gsub('_rh_', '_diff_', region_rh))
    tmp
  }
  getAveDiffOth <- function(region_rh) {
    modal <- strsplit(region_rh, '_')[[1]][1]
    region_lh <- gsub('_rh_', '_lh_', region_rh)
    region_rh_vol <- gsub(modal, 'vol', region_rh)
    region_lh_vol <- gsub('_rh_', '_lh_', region_rh_vol)
    tmp <- data.frame(
      ave=(img_df[, region_rh]*img_df[, region_rh_vol] +
        img_df[, region_lh]*img_df[, region_lh_vol])/(img_df[, region_rh_vol] + img_df[, region_lh_vol]),
      diff=(img_df[, region_rh] - img_df[, region_lh])
    )
    names(tmp) <- c(gsub('_rh_', '_ave_', region_rh), gsub('_rh_', '_diff_', region_rh))
    tmp
  }
  getLobeVol <- function(lobe) {
    vol_regions <- c(paste0('vol_rh_', regionlobe_df[regionlobe_df$lobe == lobe, 'region']),
      paste0('vol_lh_', regionlobe_df[regionlobe_df$lobe == lobe, 'region']))
    rowSums(img_df[, vol_regions])
  }
  getLobeOth <- function(lobe, modal) {
    oth_regions <- c(paste0(modal, '_rh_', regionlobe_df[regionlobe_df$lobe == lobe, 'region']),
      paste0(modal, '_lh_', regionlobe_df[regionlobe_df$lobe == lobe, 'region']))
    total_vol <- getLobeVol(lobe)
    lobe_val <- rep(0, nrow(img_df))
    for (reg in regions) {
      lobe_val <- lobe_val +
        img_df[, paste0(modal, '_rh_', reg)]*img_df[, paste0('vol_rh_', reg)]/total_vol +
        img_df[, paste0(modal, '_lh_', reg)]*img_df[, paste0('vol_lh_', reg)]/total_vol
    }
    lobe_val
  }

  vol_avediff <- do.call(bind_cols, data.frame(lapply(regions_vol_rh, getAveDiffVol)))
  ct_avediff <- do.call(bind_cols, data.frame(lapply(regions_ct_rh, getAveDiffOth)))
  gmd_avediff <- do.call(bind_cols, data.frame(lapply(regions_gmd_rh, getAveDiffOth)))

  vol_lobe <- getLobeVol(lobe)
  ct_lobe <- getLobeOth(lobe, modal='ct')
  gmd_lobe <- getLobeOth(lobe, modal='gmd')

  tmp <- cbind(vol_avediff, ct_avediff, gmd_avediff, vol_lobe, ct_lobe, gmd_lobe)
  names(tmp) <- c(names(vol_avediff), names(ct_avediff), names(gmd_avediff),
    paste0('vol_', lobe), paste0('ct_', lobe), paste0('gmd_', lobe))
  names(tmp) <- gsub('\\..*([0-9]).*', '', names(tmp))
  tmp
}

################################# Read in data #################################

img_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/imaging/antslong_struc_2021-04-26.csv')
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

plotcols <- grep('mprage', names(img_df), value=TRUE)
names(img_df) <- gsub('mprage_jlf_', '', names(img_df))
plotcols <- gsub('mprage_jlf_', '', plotcols)

region <- sapply(strsplit(plotcols, '_'), `[[`, 3)
regionlobe_df <- data.frame(region=unique(region), lobe=c('frontal', 'frontal',
  'occipital', 'temporal', 'temporal', 'parietal', 'temporal', 'parietal',
  'occipital', 'frontal', 'occipital', 'frontal', 'temporal', 'temporal',
  'frontal', 'frontal', 'frontal', 'frontal', 'occipital', 'parietal',
  'parietal', 'frontal', 'parietal', 'frontal', 'frontal', 'frontal', 'parietal',
  'temporal', 'parietal', 'temporal', 'other'))

diags <- as.character(unique(img_df$Diagnoses)[!(unique(img_df$Diagnoses) %in% c('TD-TD', 'PS-PS'))])

lobes <- unique(regionlobe_df$lobe)
derived_df <- do.call(bind_cols, sapply(lobes, getMoreRegions))

img_df <- cbind(img_df, derived_df)

plotcols <- names(derived_df)[!(names(derived_df) %in% c(grep('other', names(derived_df), value=TRUE), grep('insula', names(derived_df), value=TRUE)))]
plotcols <- plotcols[!(plotcols %in% c(grep('ave', plotcols, value=TRUE), grep('diff', plotcols, value=TRUE)))]

i=1
for (Value in plotcols) {
  img_Value_df <- img_df
  img_Value_df <- img_Value_df[!is.na(img_Value_df[, Value]), ]
  row.names(img_Value_df) <- 1:nrow(img_Value_df)

  mod1b <- gamm4(as.formula(paste(Value, "~ Diagnoses + s(Age, k=4, bs='cr') +
    s(Age, by=oDiagnoses, k=4, bs='cr')")), knots=list(15, 18, 21, 24),
    data=img_Value_df, random=~(1|bblid), REML=TRUE)

  mod2b <- gamm4(as.formula(paste(Value, "~ Male + White + Diagnoses + s(Age, k=4, bs='cr') +
    s(Age, by=oDiagnoses, k=4, bs='cr')")), knots=list(15, 18, 21, 24),
    data=img_Value_df, random=~(1|bblid), REML=TRUE)

  print(tab_model(mod1b$gam, file=paste0('~/Documents/pncLongitudinalPsychosis/results/imaging/table_', Value, '_antsknots.html')))
  assign(paste0(Value, '_model'), mod1b$gam)
  print(tab_model(mod2b$gam, file=paste0('~/Documents/pncLongitudinalPsychosis/results/imaging/tableSensitivity_', Value, '_antsknots.html')))
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

    pdf(file=paste0('~/Documents/pncLongitudinalPsychosis/plots/imaging/', Value, '_', group, '_antsknots.pdf'), width=4, height=4)
    print(img_plot)
    dev.off()
  }
  i=i+1
}








#
