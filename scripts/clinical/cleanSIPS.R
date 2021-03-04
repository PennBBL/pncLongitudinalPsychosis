### This script filters for only probands, calculates the age at each assessment,
### and sums each of the SIPS domains
###
### Ellyn Butler
### October 5, 2020 (February 12, 2021: fix two assessor problem)

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

# Identify rows with the same bblid and date
redundantRows <- function(i) {
  bblid <- final_df[i, 'bblid']
  dosips <- final_df[i, 'dosips']
  tmp_df <- final_df[final_df$bblid == bblid & final_df$dosips == dosips, ]
  if (nrow(tmp_df) > 1) {
    # Count the number of '.' in each row, keep the one with the fewest
    rows_list <- list(NA, NA)
    names(rows_list) <- c(row.names(tmp_df))
    tmp_df$numna <- rowSums(is.na(tmp_df))
    if (i == row.names(tmp_df[which.min(tmp_df$numna),])) {
      FALSE
    } else { TRUE }
    # If each rows has the same number of '.', keep the first one
  } else {
    FALSE
  }
}


# Filter data
final_df <- final_df[, c('bblid', 'dosips', 'age', 'sex', 'race', 'gaf_c', 'gaf_h',
  'first_deg_scz', 'Positive', 'Negative', 'Disorganized', 'General')]
final_df$sex <- recode(final_df$sex, `2`='Female', `1`='Male')

final_df$redundant <- sapply(1:nrow(final_df), redundantRows)
final_df <- final_df[final_df$redundant == FALSE, ]


# Export
write.csv(final_df, '~/Documents/pncLongitudinalPsychosis/data/clinical/clean_sips_nodup.csv', row.names=FALSE)
