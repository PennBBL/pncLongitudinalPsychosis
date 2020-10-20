### This script checks if there are problems with the timepoint variable
### in the raw repeated measures CNB data
###
### Ellyn Butler
### October 19, 2020

pmat_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/cognitive/cnb_practice_effects/pmat_20191113.csv')
cpw_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/cognitive/cnb_practice_effects/cpw_20191113.csv')

pmat_df[pmat_df$bblid == 80688 & pmat_df$timepoint == 1,]
cpw_df[cpw_df$bblid == 80688 & cpw_df$timepoint == 1,]

# ^ timepoint variable is wrong
