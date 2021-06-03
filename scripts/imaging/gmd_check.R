### This script correlates the XCP circa 2015 version of GMD with Ellyn's
### longitudinal pipeline's implementation. This was done because in the full
### PNC GO1 sample/XCP circa 2015, GMD increases with age, but in the longitudinal
### sample with Ellyn's method, it decreases, so we want to know if they correspond
### at all. Note that the two pipelines did not use the same parcellations.
###
### Ellyn Butler
### June 3, 2021

old_df <- read.csv('~/Documents/hiLo/data/nomeanLR/gmdData.csv')
new_df <- read.csv('~/Documents/pncLongitudinalPsychosis/data/imaging/antslong_struc_2021-04-26.csv')
demo_df <- read.csv('~/Documents/ExtraLong/data/demographicsClinical/scanid_to_seslabel_demo_20200531.csv')

new_df <- merge(demo_df, new_df)
new_df <- new_df[new_df$timepoint == 1, ]
row.names(new_df) <- 1:nrow(new_df)

df <- merge(new_df, old_df)

cor(df$mprage_jlf_gmd_R_Ent, df$mprage_jlf_gmd_rh_entorhinal) # 0.085

cor(df$mprage_jlf_gmd_R_Cun, df$mprage_jlf_vol_rh_cuneus) # -0.119

cor(df$mprage_jlf_gmd_R_FuG, df$mprage_jlf_gmd_rh_fusiform) # 0.038

cor(df$mprage_jlf_gmd_R_TTG, df$mprage_jlf_gmd_rh_transversetemporal) # -0.046

cor(df$mprage_jlf_gmd_R_PrG, df$mprage_jlf_gmd_rh_precentral) # 0.024
