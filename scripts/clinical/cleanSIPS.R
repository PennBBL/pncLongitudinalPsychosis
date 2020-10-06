### This script filters for only probands, calculates the age at each assessment,
### and sums each of the SIPS domains
###
### Ellyn Butler
### October 5, 2020

library('dplyr')
library('lubridate')

clin_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/clinical/SIPS_n752_202007.csv')
dob_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/demographics/baseline/n9498_demo_sex_race_ethnicity_dob.csv')

names(clin_df) <- tolower(names(clin_df))
clin_df <- clin_df[clin_df$type == 'Proband', ] #728 unique bblids
row.names(clin_df) <- 1:nrow(clin_df)

final_df <- merge(clin_df, dob_df, by='bblid')
final_df$dosips <- as.Date(final_df$dosips, '%m/%d/%y')
final_df$dob <- as.Date(final_df$dob, '%m/%d/%y')
final_df$age <- lubridate::time_length(difftime(final_df$dosips, final_df$dob), 'years')

pos <- paste0('p', 1:5)
neg <- paste0('n', 1:6)
dis <- paste0('d', 1:4)
gen <- paste0('g', 1:4)
gaf <- c('gaf_c', 'gaf_h')
all_domains <- c(pos, neg, dis, gen, gaf)
final_df[all_domains] <- sapply(final_df[all_domains], as.character)

# Recode . and 9 to missing
for (dom in all_domains) {
  final_df[, dom] <- na_if(final_df[, dom], '9')
  final_df[, dom] <- na_if(final_df[, dom], '.')
  final_df[, dom] <- as.numeric(final_df[, dom])
}

# Calculate sum scores
final_df$Positive <- rowSums(final_df[, pos])
final_df$Negative <- rowSums(final_df[, neg])
final_df$Disorganized <- rowSums(final_df[, dis])
final_df$General <- rowSums(final_df[, gen])

# Filter and export data
final_df <- final_df[, c('bblid', 'age', 'sex', 'race', 'gaf_c', 'gaf_h',
  'first_deg_scz', 'Positive', 'Negative', 'Disorganized', 'General')]
final_df$sex <- recode(final_df$sex, `2`='Female', `1`='Male')

write.csv(final_df, '~/Documents/pncLongitudinalPsychosis/data/clinical/clean_sips.csv', row.names=FALSE)
