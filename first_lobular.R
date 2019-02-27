### This script produces plots of lobular cortical thickness, gray matter density, gray matter volume, 
### and surface area across development, split by sex and diagnosis (second time point), and matched on race
### 
### February 15, 2019 - present

# Load libraries
library('ggplot2')
library('reshape2')
library('gridExtra')
library('lme4')
library('pbkrtest')

# Source my functions
source("/home/butellyn/ButlerPlotFuncs/plotFuncs_AdonPath.R")

# Load data
surf <- read.csv("/home/butellyn/longitudinal_psychosis/data/n2416_freesurferSurfaceArea_20170302.csv")
surf$Frontal_SA <- surf$lh_superiorfrontal_area + surf$rh_superiorfrontal_area + surf$lh_caudalmiddlefrontal_area + surf$rh_caudalmiddlefrontal_area + surf$lh_parsopercularis_area + surf$lh_parsorbitalis_area + surf$lh_parstriangularis_area + surf$rh_parsopercularis_area + surf$rh_parsorbitalis_area + surf$rh_parstriangularis_area + surf$lh_lateralorbitofrontal_area + surf$lh_medialorbitofrontal_area + surf$rh_lateralorbitofrontal_area + surf$rh_medialorbitofrontal_area + surf$lh_precentral_area + surf$rh_precentral_area + surf$lh_paracentral_area + surf$rh_paracentral_area + surf$lh_frontalpole_area + surf$rh_frontalpole_area

surf$Parietal_SA <- surf$lh_superiorparietal_area + surf$rh_superiorparietal_area + surf$lh_inferiorparietal_area + surf$rh_inferiorparietal_area + surf$lh_supramarginal_area + surf$rh_supramarginal_area + surf$lh_postcentral_area + surf$rh_postcentral_area + surf$lh_precuneus_area + surf$rh_precuneus_area

surf$Temporal_SA <- surf$lh_superiortemporal_area + surf$rh_superiortemporal_area + surf$lh_middletemporal_area + surf$rh_middletemporal_area + surf$lh_inferiortemporal_area + surf$rh_inferiortemporal_area + surf$lh_bankssts_area + surf$rh_bankssts_area + surf$lh_fusiform_area + surf$rh_fusiform_area + surf$lh_transversetemporal_area + surf$rh_transversetemporal_area + surf$lh_entorhinal_area + surf$rh_entorhinal_area + surf$lh_temporalpole_area + surf$rh_temporalpole_area + surf$lh_parahippocampal_area + surf$rh_parahippocampal_area

surf$Occipital_SA <- surf$lh_lateraloccipital_area + surf$rh_lateraloccipital_area + surf$lh_lingual_area + surf$rh_lingual_area + surf$lh_cuneus_area + surf$rh_cuneus_area + surf$lh_pericalcarine_area + surf$rh_pericalcarine_area

surf <- surf[,c("bblid", "scanid", "Frontal_SA", "Temporal_SA", "Parietal_SA", "Occipital_SA")] 



clinical <- read.csv("/home/butellyn/longitudinal_psychosis/data/pnc_diagnosis_categorical_20170526.csv")

# match on date
demo <- read.csv("/home/butellyn/longitudinal_psychosis/data/n2416_demographics_20170310.csv")
demo <- demo[,c("bblid", "scanid", "scanageMonths", "sex", "race", "DOSCAN")]
demo$DOSCAN <- as.Date(demo$DOSCAN, format = "%m/%d/%y")

clinical$bblid <- as.factor(clinical$bblid)
clinical$dodiagnosis_t1 <- as.Date(clinical$dodiagnosis_t1, format = "%m/%d/%y")
clinical$dodiagnosis_t2 <- as.Date(clinical$dodiagnosis_t2, format = "%m/%d/%y")
clinical$dodiagnosis_t3 <- as.Date(clinical$dodiagnosis_t3, format = "%m/%d/%y")
clinical$dodiagnosis_t4 <- as.Date(clinical$dodiagnosis_t4, format = "%m/%d/%y")

# PS vs. TD Terminal
demo$dx_ps <- NA
for (row in 1:nrow(demo)) {
	doscan <- demo[row, "DOSCAN"]
	bblid <- demo[row, "bblid"]
	dodiagnoses <- clinical[clinical$bblid == bblid, c("dodiagnosis_t1", "dodiagnosis_t2", "dodiagnosis_t3", "dodiagnosis_t4")]
	if (dim(dodiagnoses)[1] != 0) {
		if (!(is.na(dodiagnoses[["dodiagnosis_t1"]]))) { 
			doscan_to_dodiag_t1 <- abs(as.numeric(as.Date(doscan, format="%Y-%m-%d") - as.Date(dodiagnoses[["dodiagnosis_t1"]], format="%Y-%m-%d")))
		} else { doscan_to_dodiag_t1 <- NA }
		if (!(is.na(dodiagnoses[["dodiagnosis_t2"]]))) { 
			doscan_to_dodiag_t2 <- abs(as.numeric(as.Date(doscan, format="%Y-%m-%d") - as.Date(dodiagnoses[["dodiagnosis_t2"]], format="%Y-%m-%d")))
		} else { doscan_to_dodiag_t2 <- NA }
		if (!(is.na(dodiagnoses[["dodiagnosis_t3"]]))) { 
			doscan_to_dodiag_t3 <- abs(as.numeric(as.Date(doscan, format="%Y-%m-%d") - as.Date(dodiagnoses[["dodiagnosis_t3"]], format="%Y-%m-%d")))
		} else { doscan_to_dodiag_t3 <- NA }
		if (!(is.na(dodiagnoses[["dodiagnosis_t4"]]))) { 
			doscan_to_dodiag_t4 <- abs(as.numeric(as.Date(doscan, format="%Y-%m-%d") - as.Date(dodiagnoses[["dodiagnosis_t4"]], format="%Y-%m-%d")))
		} else { doscan_to_dodiag_t4 <- NA }
		
		# pt
		diagnosiscol <- "dx_t1_psychopathology"
		if (!(is.na(dodiagnoses[["dodiagnosis_t2"]]))) { if (doscan_to_dodiag_t1 > doscan_to_dodiag_t2) { diagnosiscol <- "dx_t2_psychopathology" }}
		if (!(is.na(dodiagnoses[["dodiagnosis_t3"]]))) { if (doscan_to_dodiag_t2 > doscan_to_dodiag_t3) { diagnosiscol <- "dx_t3_psychopathology" }}
		if (!(is.na(dodiagnoses[["dodiagnosis_t4"]]))) { if (doscan_to_dodiag_t3 > doscan_to_dodiag_t4) { diagnosiscol <- "dx_t4_psychopathology" }}
		
		diagnosis_pt <- as.character(clinical[clinical$bblid == bblid, diagnosiscol])

		# ps
		diagnosiscol <- "dx_t1_psychosis"
		if (!(is.na(dodiagnoses[["dodiagnosis_t2"]]))) { if (doscan_to_dodiag_t1 > doscan_to_dodiag_t2) { diagnosiscol <- "dx_t2_psychosis" }}
		if (!(is.na(dodiagnoses[["dodiagnosis_t3"]]))) { if (doscan_to_dodiag_t2 > doscan_to_dodiag_t3) { diagnosiscol <- "dx_t3_psychosis" }}
		if (!(is.na(dodiagnoses[["dodiagnosis_t4"]]))) { if (doscan_to_dodiag_t3 > doscan_to_dodiag_t4) { diagnosiscol <- "dx_t4_psychosis" }}
		
		diagnosis_ps <- as.character(clinical[clinical$bblid == bblid, diagnosiscol])

		# change PTs to PSs where appropriate ####### KOSHA SHOULD I DO THIS?
		if (diagnosis_pt == "PT" & diagnosis_ps == "PS") { diagnosis <- "PS" 
		} else if (diagnosis_pt == "PT" & diagnosis_ps != "PS") { diagnosis <- "OP"
		} else if (diagnosis_pt == "PT" & diagnosis_ps == "NONPS") { print(paste(bblid, "eek!"))
		} else if (diagnosis_pt == "TD" & diagnosis_ps == "NONPS") { diagnosis <- "TD" }

		demo[row, "dx_ps"] <- diagnosis
	}
}
demo$dx_ps <- as.factor(demo$dx_ps)

demo$Age <- demo$scanageMonths/12
demo$bblid <- as.factor(demo$bblid)
for (bblid in levels(demo$bblid)) {
	tmp <- demo[demo$bblid == bblid, ]
	lasttimerow <- rownames(tmp[tmp$Age == max(tmp$Age),])
	lastdiag <- as.character(tmp[lasttimerow, "dx_ps"])
	demo[demo$bblid == bblid, "psTerminal"] <- lastdiag
}
demo$psTerminal <- as.factor(demo$psTerminal)

# February 27, 2019: Get rid of subjects who have OP as their final diagnosis
toBeRemoved <- which(demo$psTerminal == "OP")
demo <- demo[-toBeRemoved,]
demo$psTerminal <- factor(demo$psTerminal)

# Persister vs. TD # February 26, 2019: Something wrong here 
#demo$ps_trajectory <- NA
#for (row in 1:nrow(demo)) {
#	bblid <- demo[row, "bblid"]
#	ps_trajectory <- as.character(clinical[clinical$bblid == bblid, "pncGrpPsychosisCl"])
#	if (length(ps_trajectory) != 0) {
#		if (ps_trajectory == "TD") {
#			demo[row, "ps_trajectory"] <- ps_trajectory
#		} else if (ps_trajectory == "Persister") {
#			demo[row, "ps_trajectory"] <- ps_trajectory
#		}
#	}
#}
#demo$ps_trajectory <- as.factor(demo$ps_trajectory)

# start imaging data
cort <- read.csv("/home/butellyn/longitudinal_psychosis/data/n2416_jlfAntsCTIntersectionCt_20170331.csv")
vol <- read.csv("/home/butellyn/longitudinal_psychosis/data/n2416_jlfAntsCTIntersectionVol_20170412.csv")
gmd <- read.csv("/home/butellyn/longitudinal_psychosis/data/n2416_jlfAtroposIntersectionGMD_20170410.csv")

# create dfs for lobular values
volcortcols <- c("mprage_jlf_vol_3rd_Ventricle", "mprage_jlf_vol_4th_Ventricle", "mprage_jlf_vol_R_Accumbens_Area", "mprage_jlf_vol_L_Accumbens_Area", "mprage_jlf_vol_R_Amygdala", "mprage_jlf_vol_L_Amygdala", "mprage_jlf_vol_Brain_Stem", "mprage_jlf_vol_R_Caudate", "mprage_jlf_vol_L_Caudate", "mprage_jlf_vol_R_Cerebellum_Exterior", "mprage_jlf_vol_L_Cerebellum_Exterior", "mprage_jlf_vol_R_Cerebellum_White_Matter", "mprage_jlf_vol_L_Cerebellum_White_Matter", "mprage_jlf_vol_R_Cerebral_White_Matter", "mprage_jlf_vol_L_Cerebral_White_Matter", "mprage_jlf_vol_CSF", "mprage_jlf_vol_R_Hippocampus", "mprage_jlf_vol_L_Hippocampus", "mprage_jlf_vol_R_Inf_Lat_Vent", "mprage_jlf_vol_L_Inf_Lat_Vent", "mprage_jlf_vol_R_Lateral_Ventricle", "mprage_jlf_vol_L_Lateral_Ventricle", "mprage_jlf_vol_R_Pallidum", "mprage_jlf_vol_L_Pallidum", "mprage_jlf_vol_R_Putamen", "mprage_jlf_vol_L_Putamen", "mprage_jlf_vol_R_Thalamus_Proper", "mprage_jlf_vol_L_Thalamus_Proper", "mprage_jlf_vol_Cerebellar_Vermal_Lobules_I.V", "mprage_jlf_vol_Cerebellar_Vermal_Lobules_VI.VII", "mprage_jlf_vol_Cerebellar_Vermal_Lobules_VIII.X", "mprage_jlf_vol_R_ACgG", "mprage_jlf_vol_L_ACgG", "mprage_jlf_vol_R_AIns", "mprage_jlf_vol_L_AIns", "mprage_jlf_vol_R_PHG", "mprage_jlf_vol_L_PHG", "mprage_jlf_vol_R_PIns", "mprage_jlf_vol_L_PIns", "mprage_jlf_vol_R_SCA", "mprage_jlf_vol_L_SCA", "mprage_jlf_vol_R_PCgG", "mprage_jlf_vol_L_PCgG", "mprage_jlf_vol_R_Ent", "mprage_jlf_vol_L_Ent", "mprage_jlf_vol_R_MCgG", "mprage_jlf_vol_L_MCgG")
vol_for_cort_df <- vol[,!(names(vol) %in% volcortcols)]
cortcols <- c("mprage_jlf_ct_R_PHG", "mprage_jlf_ct_L_PHG", "mprage_jlf_ct_R_PIns", "mprage_jlf_ct_L_PIns", "mprage_jlf_ct_R_SCA", "mprage_jlf_ct_L_SCA", "mprage_jlf_ct_R_AIns", "mprage_jlf_ct_L_AIns", "mprage_jlf_ct_R_ACgG", "mprage_jlf_ct_L_ACgG", "mprage_jlf_ct_R_PCgG", "mprage_jlf_ct_L_PCgG", "mprage_jlf_ct_R_Ent", "mprage_jlf_ct_L_Ent", "mprage_jlf_ct_R_MCgG", "mprage_jlf_ct_L_MCgG")
cort_df <- cort[,!(names(cort) %in% cortcols)]

volgmdcols <- c("mprage_jlf_vol_3rd_Ventricle", "mprage_jlf_vol_4th_Ventricle", "mprage_jlf_vol_R_Cerebellum_Exterior", "mprage_jlf_vol_L_Cerebellum_Exterior", "mprage_jlf_vol_R_Cerebellum_White_Matter", "mprage_jlf_vol_L_Cerebellum_White_Matter", "mprage_jlf_vol_R_Cerebral_White_Matter", "mprage_jlf_vol_L_Cerebral_White_Matter", "mprage_jlf_vol_CSF", "mprage_jlf_vol_R_Inf_Lat_Vent", "mprage_jlf_vol_L_Inf_Lat_Vent", "mprage_jlf_vol_R_Lateral_Ventricle", "mprage_jlf_vol_L_Lateral_Ventricle", "mprage_jlf_vol_Cerebellar_Vermal_Lobules_I.V", "mprage_jlf_vol_Cerebellar_Vermal_Lobules_VI.VII", "mprage_jlf_vol_Cerebellar_Vermal_Lobules_VIII.X")
vol_for_gmd_df <- vol[,!(names(vol) %in% volgmdcols)]
gmdcols <- c("mprage_jlf_gmd_Brain_Stem", "mprage_jlf_gmd_R_Cerebellum_Exterior", "mprage_jlf_gmd_L_Cerebellum_Exterior", "mprage_jlf_gmd_Cerebellar_Vermal_Lobules_I.V", "mprage_jlf_gmd_Cerebellar_Vermal_Lobules_VI.VII", "mprage_jlf_gmd_Cerebellar_Vermal_Lobules_VIII.X")
gmd_df <- gmd[,!(names(gmd) %in% gmdcols)]

# create lobular values
vol_df2 <- averageLeftAndRight_Vol(vol, rightString="_R_", leftString="_L_", averageString="_ave_")

##### cort #####
cort_df2 <- averageLeftAndRight_WeightByVol(vol_for_cort_df, cort_df, volString="_vol", otherString="_ct")
ROIlist <- grep("_ave_", colnames(cort_df2), value=TRUE)
ROIlist <- addUnderScore(ROIlist)
ROI_ListofLists <- roiLobes(ROIlist, lobeDef=FALSE)
ROI_ListofLists <- removeUnderScore(ROI_ListofLists)
ROI_ListofLists$BasGang <- NULL
ROI_ListofLists$Limbic <- NULL
ROI_ListofLists$Cerebellum <- NULL

vol_for_cort <- vol_df2
ROIlist <- grep("_ave_", colnames(vol_for_cort), value=TRUE)
ROIlist <- addUnderScore(ROIlist)
ROI_ListofLists_Vol <- roiLobes(ROIlist, lobeDef=FALSE)
ROI_ListofLists_Vol <- removeUnderScore(ROI_ListofLists_Vol)
ROI_ListofLists_Vol$BasGang <- NULL
ROI_ListofLists_Vol$Limbic <- NULL
ROI_ListofLists_Vol$Cerebellum <- NULL
vol_for_cort <- lobeVolumes(vol_for_cort, ROI_ListofLists_Vol, lastList=TRUE) 

cort_df2 <- averageROIsInLobes_WeightByVol(vol_for_cort, cort_df2, ROI_ListofLists_Vol, ROI_ListofLists, lastList=TRUE, type="cort")
cort_df3 <- cort_df2[,c("bblid", "scanid", "FrontOrb_CT" , "FrontDors_CT", "Temporal_CT", "Parietal_CT", "Occipital_CT")]

# gmd
gmd_df2 <- averageLeftAndRight_WeightByVol(vol_for_gmd_df, gmd_df, volString="_vol", otherString="_gmd")
ROIlist <- grep("_ave_", colnames(gmd_df2), value=TRUE)
ROIlist <- addUnderScore(ROIlist)
ROI_ListofLists <- roiLobes(ROIlist, lobeDef=FALSE)
ROI_ListofLists <- removeUnderScore(ROI_ListofLists)
ROI_ListofLists$Cerebellum <- NULL

vol_for_gmd <- vol_df2
ROIlist <- grep("_ave_", colnames(vol_for_gmd), value=TRUE)
ROIlist <- addUnderScore(ROIlist)
ROI_ListofLists_Vol <- roiLobes(ROIlist, lobeDef=FALSE)
ROI_ListofLists_Vol <- removeUnderScore(ROI_ListofLists_Vol)
ROI_ListofLists_Vol$Cerebellum <- NULL
vol_for_gmd <- lobeVolumes(vol_for_gmd, ROI_ListofLists_Vol, lastList=TRUE) 

gmd_df2 <- averageROIsInLobes_WeightByVol(vol_for_gmd, gmd_df2, ROI_ListofLists_Vol, ROI_ListofLists, lastList=TRUE, type="gmd")
gmd_df3 <- gmd_df2[,c("bblid", "scanid", "BasGang_GMD", "Limbic_GMD", "FrontOrb_GMD", "FrontDors_GMD", "Temporal_GMD", "Parietal_GMD", "Occipital_GMD")]

# vol
ROIlist <- grep("_R_", colnames(vol_df2), value=TRUE)
ROIlist <- c(ROIlist, grep("_L_", colnames(vol_df2), value=TRUE))
ROIlist <- addUnderScore(ROIlist)
ROI_ListofLists <- roiLobes(ROIlist, lobeDef=FALSE)
ROI_ListofLists <- removeUnderScore(ROI_ListofLists)

vol_df2 <- lobeVolumes(vol_df2, ROI_ListofLists)
vol_df3 <- vol_df2[,c("bblid", "scanid", "BasGang_Vol", "Limbic_Vol", "FrontOrb_Vol", "FrontDors_Vol", "Temporal_Vol", "Parietal_Vol", "Occipital_Vol")]

# bblid and scanid as factor
demo$bblid <- as.factor(demo$bblid)
surf$bblid <- as.factor(surf$bblid)
cort_df3$bblid <- as.factor(cort_df3$bblid)
gmd_df3$bblid <- as.factor(gmd_df3$bblid)
vol_df3$bblid <- as.factor(vol_df3$bblid)
demo$scanid <- as.factor(demo$scanid)
surf$scanid <- as.factor(surf$scanid)
cort_df3$scanid <- as.factor(cort_df3$scanid)
gmd_df3$scanid <- as.factor(gmd_df3$scanid)
vol_df3$scanid <- as.factor(vol_df3$scanid)

df <- merge(demo, surf, by=c("bblid","scanid"))
df <- merge(df, cort_df3, by=c("bblid","scanid"))
df <- merge(df, gmd_df3, by=c("bblid","scanid"))
df <- merge(df, vol_df3, by=c("bblid","scanid"))


df$sex <- as.factor(df$sex)
df$Age <- df$scanageMonths/12

#### PS
# get rid of rows without diagnosis
df <- df[!(is.na(df$psTerminal)),]

df$bblid <- factor(df$bblid)

# Get rid of people with only one time point
delrows <- c()
for (bblid in levels(df$bblid)) {
	tmp <- df[df$bblid == bblid, ]
	if (dim(tmp)[1] == 1) {
		delrows <- c(delrows, rownames(tmp))
	}
}
df <- df[!(rownames(df) %in% delrows),]


####################################### Start analyses #######################################

N_F_ps <- nrow(df[df$sex == 2 & df$psTerminal == "PS" & !(duplicated(df["bblid"])),])
per_W_F_ps <- 100*round(nrow(df[df$sex == 2 & df$psTerminal == "PS" & df$race == 1 & !(duplicated(df["bblid"])),])/N_F_ps, digits=3)
per_B_F_ps <- 100*round(nrow(df[df$sex == 2 & df$psTerminal == "PS" & df$race == 2 & !(duplicated(df["bblid"])),])/N_F_ps, digits=3)

N_F_td <- nrow(df[df$sex == 2 & df$psTerminal == "TD" & !(duplicated(df["bblid"])),])
per_W_F_td <- 100*round(nrow(df[df$sex == 2 & df$psTerminal == "TD" & df$race == 1 & !(duplicated(df["bblid"])),])/N_F_td, digits=3)
per_B_F_td <- 100*round(nrow(df[df$sex == 2 & df$psTerminal == "TD" & df$race == 2 & !(duplicated(df["bblid"])),])/N_F_td, digits=3)

#N_F_pt <- nrow(df[df$sex == 2 & df$psTerminal == "OP" & !(duplicated(df["bblid"])),])
#per_W_F_pt <- 100*round(nrow(df[df$sex == 2 & df$psTerminal == "OP" & df$race == 1 & !(duplicated(df["bblid"])),])/N_F_pt, digits=3)
#per_B_F_pt <- 100*round(nrow(df[df$sex == 2 & df$psTerminal == "OP" & df$race == 2 & !(duplicated(df["bblid"])),])/N_F_pt, digits=3)

N_M_ps <- nrow(df[df$sex == 1 & df$psTerminal == "PS" & !(duplicated(df["bblid"])),])
per_W_M_ps <- 100*round(nrow(df[df$sex == 1 & df$psTerminal == "PS" & df$race == 1 & !(duplicated(df["bblid"])), ])/N_M_ps, digits=3)
per_B_M_ps <- 100*round(nrow(df[df$sex == 1 & df$psTerminal == "PS" & df$race == 2 & !(duplicated(df["bblid"])),])/N_M_ps, digits=3)

N_M_td <- nrow(df[df$sex == 1 & df$psTerminal == "TD" & !(duplicated(df["bblid"])),])
per_W_M_td <- 100*round(nrow(df[df$sex == 1 & df$psTerminal == "TD" & df$race == 1 & !(duplicated(df["bblid"])),])/N_M_td, digits=3)
per_B_M_td <- 100*round(nrow(df[df$sex == 1 & df$psTerminal == "TD" & df$race == 2 & !(duplicated(df["bblid"])),])/N_M_td, digits=3)

#N_M_pt <- nrow(df[df$sex == 1 & df$psTerminal == "OP" & !(duplicated(df["bblid"])),])
#per_W_M_pt <- 100*round(nrow(df[df$sex == 1 & df$psTerminal == "OP" & df$race == 1 & !(duplicated(df["bblid"])),])/N_M_pt, digits=3)
#per_B_M_pt <- 100*round(nrow(df[df$sex == 1 & df$psTerminal == "OP" & df$race == 2 & !(duplicated(df["bblid"])),])/N_M_pt, digits=3)

femalestring <- paste0("PS (N=", N_F_ps, ", W=", per_W_F_ps, "%, B=", per_B_F_ps, "%), TD (N=", N_F_td, ", W=", per_W_F_td, "%, B=", per_B_F_td, "%)")
malestring <- paste0("PS (N=", N_M_ps, ", W=", per_W_M_ps, "%, B=", per_B_M_ps, "%), TD (N=", N_M_td, ", W=", per_W_M_td, "%, B=", per_B_M_td, "%)")


#----------------------------------- Female -----------------------------------#
# Frontal Plots
stats_frontorb_cort_F <- lmer(FrontOrb_CT ~ Age*psTerminal + (Age | bblid), df[df$sex == 2,])
#smallMod_noAge <- lmer(FrontOrb_CT ~ psTerminal + Age:psTerminal + (Age | bblid), df[df$sex == 2,]) # February 27, 2019: getting two interaction terms
#sig_age <- KRmodcomp(stats_frontorb_cort_F, smallMod_noAge, betaH = 0, details = 0)
#smallMod_noDiag <- lmer(FrontOrb_CT ~ Age + Age:psTerminal + (Age | bblid), df[df$sex == 2,])
#smallMod_noInter <- lmer(FrontOrb_CT ~ Age + psTerminal + (Age | bblid), df[df$sex == 2,])
coeffdf <- summary(stats_frontorb_cort_F)$coefficients
estimates <- coeffdf[,"Estimate"]
#stderrors <- coeffdf[,"Std. Error"]
tvalues <- coeffdf[,"t value"]
if (abs(tvalues[["Age"]]) > 2) { age_est <- paste0(round(estimates[["Age"]], digits=3), "*") } else { age_est <- round(estimates[["Age"]], digits=3) }
if (abs(tvalues[["psTerminalTD"]]) > 2) { diag_est <- paste0(round(estimates[["psTerminalTD"]], digits=3), "*") } else { diag_est <- round(estimates[["psTerminalTD"]], digits=3) }
if (abs(tvalues[["Age:psTerminalTD"]]) > 2) { inter_est <- paste0(round(estimates[["Age:psTerminalTD"]], digits=3), "*") } else { inter_est <- round(estimates[["Age:psTerminalTD"]], digits=3) }
mainstring <- paste0(femalestring, "\n", "Coefficients: Age=", age_est, ", Diagnosis=", diag_est, ", Age:Diagnosis=", inter_est)
frontorb_cort_F <- ggplot(df[df$sex == 2,], aes(x=Age, y=FrontOrb_CT)) + 
	geom_line(aes(group=bblid, col=psTerminal), alpha=.2) + 
	geom_smooth(se = TRUE, method = "glm", aes(group=psTerminal, col=psTerminal)) + 
	theme_minimal() + scale_colour_manual(values = c("mediumaquamarine", "maroon4")) +
	theme(plot.title = element_text(family="Times", face="bold", size=20)) +
	ylab("Cortical Thickness (mm)") +
	labs(title = "Frontal Orbital Cortical Thickness (Female)", subtitle = mainstring)

stats_frontorb_vol_F <- lmer(FrontOrb_Vol ~ Age*psTerminal + (Age | bblid), df[df$sex == 2,])
coeffdf <- summary(stats_frontorb_vol_F)$coefficients
estimates <- coeffdf[,"Estimate"]
stderrors <- coeffdf[,"Std. Error"]
tvalues <- coeffdf[,"t value"]
if (abs(tvalues[["Age"]]) > 2) { age_est <- paste0(round(estimates[["Age"]], digits=3), "*") } else { age_est <- round(estimates[["Age"]], digits=3) }
if (abs(tvalues[["psTerminalTD"]]) > 2) { diag_est <- paste0(round(estimates[["psTerminalTD"]], digits=3), "*") } else { diag_est <- round(estimates[["psTerminalTD"]], digits=3) }
if (abs(tvalues[["Age:psTerminalTD"]]) > 2) { inter_est <- paste0(round(estimates[["Age:psTerminalTD"]], digits=3), "*") } else { inter_est <- round(estimates[["Age:psTerminalTD"]], digits=3) }
mainstring <- paste0(femalestring, "\n", "Coefficients: Age=", age_est, ", Diagnosis=", diag_est, ", Age:Diagnosis=", inter_est)
frontorb_vol_F <- ggplot(df[df$sex == 2,], aes(x=Age, y=FrontOrb_Vol)) + 
	geom_line(aes(group=bblid, col=psTerminal), alpha=.2) + 
	geom_smooth(se = TRUE, method = "glm", aes(group=psTerminal, col=psTerminal)) + 
	theme_minimal() + scale_colour_manual(values = c("mediumaquamarine", "maroon4")) +
	theme(plot.title = element_text(family="Times", face="bold", size=20)) +
	ylab("Gray Matter Volume (mm3)") +
	labs(title = "Frontal Orbital Gray Matter Volume (Female)", subtitle = mainstring)

stats_frontorb_gmd_F <- lmer(FrontOrb_GMD ~ Age*psTerminal + (Age | bblid), df[df$sex == 2,])
coeffdf <- summary(stats_frontorb_gmd_F)$coefficients
estimates <- coeffdf[,"Estimate"]
stderrors <- coeffdf[,"Std. Error"]
tvalues <- coeffdf[,"t value"]
if (abs(tvalues[["Age"]]) > 2) { age_est <- paste0(round(estimates[["Age"]], digits=3), "*") } else { age_est <- round(estimates[["Age"]], digits=3) }
if (abs(tvalues[["psTerminalTD"]]) > 2) { diag_est <- paste0(round(estimates[["psTerminalTD"]], digits=3), "*") } else { diag_est <- round(estimates[["psTerminalTD"]], digits=3) }
if (abs(tvalues[["Age:psTerminalTD"]]) > 2) { inter_est <- paste0(round(estimates[["Age:psTerminalTD"]], digits=3), "*") } else { inter_est <- round(estimates[["Age:psTerminalTD"]], digits=3) }
mainstring <- paste0(femalestring, "\n", "Coefficients: Age=", age_est, ", Diagnosis=", diag_est, ", Age:Diagnosis=", inter_est)
frontorb_gmd_F <- ggplot(df[df$sex == 2,], aes(x=Age, y=FrontOrb_GMD)) + 
	geom_line(aes(group=bblid, col=psTerminal), alpha=.2) + 
	geom_smooth(se = TRUE, method = "glm", aes(group=psTerminal, col=psTerminal)) + 
	theme_minimal() + scale_colour_manual(values = c("mediumaquamarine", "maroon4")) +
	theme(plot.title = element_text(family="Times", face="bold", size=20)) +
	ylab("Gray Matter Density") +
	labs(title = "Frontal Orbital Gray Matter Density (Female)", subtitle = mainstring)

stats_frontdors_cort_F <- lmer(FrontDors_CT ~ Age*psTerminal + (Age | bblid), df[df$sex == 2,])
coeffdf <- summary(stats_frontdors_cort_F)$coefficients
estimates <- coeffdf[,"Estimate"]
stderrors <- coeffdf[,"Std. Error"]
tvalues <- coeffdf[,"t value"]
if (abs(tvalues[["Age"]]) > 2) { age_est <- paste0(round(estimates[["Age"]], digits=3), "*") } else { age_est <- round(estimates[["Age"]], digits=3) }
if (abs(tvalues[["psTerminalTD"]]) > 2) { diag_est <- paste0(round(estimates[["psTerminalTD"]], digits=3), "*") } else { diag_est <- round(estimates[["psTerminalTD"]], digits=3) }
if (abs(tvalues[["Age:psTerminalTD"]]) > 2) { inter_est <- paste0(round(estimates[["Age:psTerminalTD"]], digits=3), "*") } else { inter_est <- round(estimates[["Age:psTerminalTD"]], digits=3) }
mainstring <- paste0(femalestring, "\n", "Coefficients: Age=", age_est, ", Diagnosis=", diag_est, ", Age:Diagnosis=", inter_est)
frontdors_cort_F <- ggplot(df[df$sex == 2,], aes(x=Age, y=FrontDors_CT)) + 
	geom_line(aes(group=bblid, col=psTerminal), alpha=.2) + 
	geom_smooth(se = TRUE, method = "glm", aes(group=psTerminal, col=psTerminal)) + 
	theme_minimal() + scale_colour_manual(values = c("mediumaquamarine", "maroon4")) +
	theme(plot.title = element_text(family="Times", face="bold", size=20)) +
	ylab("Cortical Thickness (mm)") +
	labs(title = "Frontal Dorsal Cortical Thickness (Female)", subtitle = mainstring)

stats_frontdors_vol_F <- lmer(FrontDors_Vol ~ Age*psTerminal + (Age | bblid), df[df$sex == 2,])
coeffdf <- summary(stats_frontdors_vol_F)$coefficients
estimates <- coeffdf[,"Estimate"]
stderrors <- coeffdf[,"Std. Error"]
tvalues <- coeffdf[,"t value"]
if (abs(tvalues[["Age"]]) > 2) { age_est <- paste0(round(estimates[["Age"]], digits=3), "*") } else { age_est <- round(estimates[["Age"]], digits=3) }
if (abs(tvalues[["psTerminalTD"]]) > 2) { diag_est <- paste0(round(estimates[["psTerminalTD"]], digits=3), "*") } else { diag_est <- round(estimates[["psTerminalTD"]], digits=3) }
if (abs(tvalues[["Age:psTerminalTD"]]) > 2) { inter_est <- paste0(round(estimates[["Age:psTerminalTD"]], digits=3), "*") } else { inter_est <- round(estimates[["Age:psTerminalTD"]], digits=3) }
mainstring <- paste0(femalestring, "\n", "Coefficients: Age=", age_est, ", Diagnosis=", diag_est, ", Age:Diagnosis=", inter_est)
frontdors_vol_F <- ggplot(df[df$sex == 2,], aes(x=Age, y=FrontDors_Vol)) + 
	geom_line(aes(group=bblid, col=psTerminal), alpha=.2) + 
	geom_smooth(se = TRUE, method = "glm", aes(group=psTerminal, col=psTerminal)) + 
	theme_minimal() + scale_colour_manual(values = c("mediumaquamarine", "maroon4")) +
	theme(plot.title = element_text(family="Times", face="bold", size=20)) +
	ylab("Gray Matter Volume (mm3)") +
	labs(title = "Frontal Dorsal Gray Matter Volume (Female)", subtitle = mainstring)

stats_frontdors_gmd_F <- lmer(FrontDors_GMD ~ Age*psTerminal + (Age | bblid), df[df$sex == 2,])
coeffdf <- summary(stats_frontdors_gmd_F)$coefficients
estimates <- coeffdf[,"Estimate"]
stderrors <- coeffdf[,"Std. Error"]
tvalues <- coeffdf[,"t value"]
if (abs(tvalues[["Age"]]) > 2) { age_est <- paste0(round(estimates[["Age"]], digits=3), "*") } else { age_est <- round(estimates[["Age"]], digits=3) }
if (abs(tvalues[["psTerminalTD"]]) > 2) { diag_est <- paste0(round(estimates[["psTerminalTD"]], digits=3), "*") } else { diag_est <- round(estimates[["psTerminalTD"]], digits=3) }
if (abs(tvalues[["Age:psTerminalTD"]]) > 2) { inter_est <- paste0(round(estimates[["Age:psTerminalTD"]], digits=3), "*") } else { inter_est <- round(estimates[["Age:psTerminalTD"]], digits=3) }
mainstring <- paste0(femalestring, "\n", "Coefficients: Age=", age_est, ", Diagnosis=", diag_est, ", Age:Diagnosis=", inter_est)
frontdors_gmd_F <- ggplot(df[df$sex == 2,], aes(x=Age, y=FrontDors_GMD)) + 
	geom_line(aes(group=bblid, col=psTerminal), alpha=.2) + 
	geom_smooth(se = TRUE, method = "glm", aes(group=psTerminal, col=psTerminal)) + 
	theme_minimal() + scale_colour_manual(values = c("mediumaquamarine", "maroon4")) +
	theme(plot.title = element_text(family="Times", face="bold", size=20)) +
	ylab("Gray Matter Density") +
	labs(title = "Frontal Dorsal Gray Matter Density (Female)", subtitle = mainstring)

stats_front_sa_F <- lmer(Frontal_SA ~ Age*psTerminal + (Age | bblid), df[df$sex == 2,])
coeffdf <- summary(stats_front_sa_F)$coefficients
estimates <- coeffdf[,"Estimate"]
stderrors <- coeffdf[,"Std. Error"]
tvalues <- coeffdf[,"t value"]
if (abs(tvalues[["Age"]]) > 2) { age_est <- paste0(round(estimates[["Age"]], digits=3), "*") } else { age_est <- round(estimates[["Age"]], digits=3) }
if (abs(tvalues[["psTerminalTD"]]) > 2) { diag_est <- paste0(round(estimates[["psTerminalTD"]], digits=3), "*") } else { diag_est <- round(estimates[["psTerminalTD"]], digits=3) }
if (abs(tvalues[["Age:psTerminalTD"]]) > 2) { inter_est <- paste0(round(estimates[["Age:psTerminalTD"]], digits=3), "*") } else { inter_est <- round(estimates[["Age:psTerminalTD"]], digits=3) }
mainstring <- paste0(femalestring, "\n", "Coefficients: Age=", age_est, ", Diagnosis=", diag_est, ", Age:Diagnosis=", inter_est)
front_sa_F <- ggplot(df[df$sex == 2,], aes(x=Age, y=Frontal_SA)) + 
	geom_line(aes(group=bblid, col=psTerminal), alpha=.2) + 
	geom_smooth(se = TRUE, method = "glm", aes(group=psTerminal, col=psTerminal)) + 
	theme_minimal() + scale_colour_manual(values = c("mediumaquamarine", "maroon4")) +
	theme(plot.title = element_text(family="Times", face="bold", size=20)) +
	ylab("Surface Area (mm2)") +
	labs(title = "Frontal Surface Area (Female)", subtitle = mainstring)

# Parietal
stats_parietal_cort_F <- lmer(Parietal_CT ~ Age*psTerminal + (Age | bblid), df[df$sex == 2,])
coeffdf <- summary(stats_parietal_cort_F)$coefficients
estimates <- coeffdf[,"Estimate"]
stderrors <- coeffdf[,"Std. Error"]
tvalues <- coeffdf[,"t value"]
if (abs(tvalues[["Age"]]) > 2) { age_est <- paste0(round(estimates[["Age"]], digits=3), "*") } else { age_est <- round(estimates[["Age"]], digits=3) }
if (abs(tvalues[["psTerminalTD"]]) > 2) { diag_est <- paste0(round(estimates[["psTerminalTD"]], digits=3), "*") } else { diag_est <- round(estimates[["psTerminalTD"]], digits=3) }
if (abs(tvalues[["Age:psTerminalTD"]]) > 2) { inter_est <- paste0(round(estimates[["Age:psTerminalTD"]], digits=3), "*") } else { inter_est <- round(estimates[["Age:psTerminalTD"]], digits=3) }
mainstring <- paste0(femalestring, "\n", "Coefficients: Age=", age_est, ", Diagnosis=", diag_est, ", Age:Diagnosis=", inter_est)
parietal_cort_F <- ggplot(df[df$sex == 2,], aes(x=Age, y=Parietal_CT)) + 
	geom_line(aes(group=bblid, col=psTerminal), alpha=.2) + 
	geom_smooth(se = TRUE, method = "glm",  aes(group=psTerminal, col=psTerminal)) + 
	theme_minimal() + scale_colour_manual(values = c("mediumaquamarine", "maroon4")) +
	theme(plot.title = element_text(family="Times", face="bold", size=20)) +
	ylab("Cortical Thickness (mm)") +
	labs(title = "Parietal Cortical Thickness (Female)", subtitle = mainstring)

stats_parietal_vol_F <- lmer(Parietal_Vol ~ Age*psTerminal + (Age | bblid), df[df$sex == 2,])
coeffdf <- summary(stats_parietal_vol_F)$coefficients
estimates <- coeffdf[,"Estimate"]
stderrors <- coeffdf[,"Std. Error"]
tvalues <- coeffdf[,"t value"]
if (abs(tvalues[["Age"]]) > 2) { age_est <- paste0(round(estimates[["Age"]], digits=3), "*") } else { age_est <- round(estimates[["Age"]], digits=3) }
if (abs(tvalues[["psTerminalTD"]]) > 2) { diag_est <- paste0(round(estimates[["psTerminalTD"]], digits=3), "*") } else { diag_est <- round(estimates[["psTerminalTD"]], digits=3) }
if (abs(tvalues[["Age:psTerminalTD"]]) > 2) { inter_est <- paste0(round(estimates[["Age:psTerminalTD"]], digits=3), "*") } else { inter_est <- round(estimates[["Age:psTerminalTD"]], digits=3) }
mainstring <- paste0(femalestring, "\n", "Coefficients: Age=", age_est, ", Diagnosis=", diag_est, ", Age:Diagnosis=", inter_est)
parietal_vol_F <- ggplot(df[df$sex == 2,], aes(x=Age, y=Parietal_Vol)) + 
	geom_line(aes(group=bblid, col=psTerminal), alpha=.2) + 
	geom_smooth(se = TRUE, method = "glm", aes(group=psTerminal, col=psTerminal)) + 
	theme_minimal() + scale_colour_manual(values = c("mediumaquamarine", "maroon4")) +
	theme(plot.title = element_text(family="Times", face="bold", size=20)) +
	ylab("Gray Matter Volume (mm3)") +
	labs(title = "Parietal Gray Matter Volume (Female)", subtitle = mainstring)

stats_parietal_gmd_F <- lmer(Parietal_GMD ~ Age*psTerminal + (Age | bblid), df[df$sex == 2,])
coeffdf <- summary(stats_parietal_gmd_F)$coefficients
estimates <- coeffdf[,"Estimate"]
stderrors <- coeffdf[,"Std. Error"]
tvalues <- coeffdf[,"t value"]
if (abs(tvalues[["Age"]]) > 2) { age_est <- paste0(round(estimates[["Age"]], digits=3), "*") } else { age_est <- round(estimates[["Age"]], digits=3) }
if (abs(tvalues[["psTerminalTD"]]) > 2) { diag_est <- paste0(round(estimates[["psTerminalTD"]], digits=3), "*") } else { diag_est <- round(estimates[["psTerminalTD"]], digits=3) }
if (abs(tvalues[["Age:psTerminalTD"]]) > 2) { inter_est <- paste0(round(estimates[["Age:psTerminalTD"]], digits=3), "*") } else { inter_est <- round(estimates[["Age:psTerminalTD"]], digits=3) }
mainstring <- paste0(femalestring, "\n", "Coefficients: Age=", age_est, ", Diagnosis=", diag_est, ", Age:Diagnosis=", inter_est)
parietal_gmd_F <- ggplot(df[df$sex == 2,], aes(x=Age, y=Parietal_GMD)) + 
	geom_line(aes(group=bblid, col=psTerminal), alpha=.2) + 
	geom_smooth(se = TRUE, method = "glm", aes(group=psTerminal, col=psTerminal)) + 
	theme_minimal() + scale_colour_manual(values = c("mediumaquamarine", "maroon4")) +
	theme(plot.title = element_text(family="Times", face="bold", size=20)) +
	ylab("Gray Matter Density") +
	labs(title = "Parietal Gray Matter Density (Female)", subtitle = mainstring)

stats_parietal_sa_F <- lmer(Parietal_SA ~ Age*psTerminal + (Age | bblid), df[df$sex == 2,])
coeffdf <- summary(stats_parietal_sa_F)$coefficients
estimates <- coeffdf[,"Estimate"]
stderrors <- coeffdf[,"Std. Error"]
tvalues <- coeffdf[,"t value"]
if (abs(tvalues[["Age"]]) > 2) { age_est <- paste0(round(estimates[["Age"]], digits=3), "*") } else { age_est <- round(estimates[["Age"]], digits=3) }
if (abs(tvalues[["psTerminalTD"]]) > 2) { diag_est <- paste0(round(estimates[["psTerminalTD"]], digits=3), "*") } else { diag_est <- round(estimates[["psTerminalTD"]], digits=3) }
if (abs(tvalues[["Age:psTerminalTD"]]) > 2) { inter_est <- paste0(round(estimates[["Age:psTerminalTD"]], digits=3), "*") } else { inter_est <- round(estimates[["Age:psTerminalTD"]], digits=3) }
mainstring <- paste0(femalestring, "\n", "Coefficients: Age=", age_est, ", Diagnosis=", diag_est, ", Age:Diagnosis=", inter_est)
parietal_sa_F <- ggplot(df[df$sex == 2,], aes(x=Age, y=Parietal_SA)) + 
	geom_line(aes(group=bblid, col=psTerminal), alpha=.2) + 
	geom_smooth(se = TRUE, method = "glm", aes(group=psTerminal, col=psTerminal)) + 
	theme_minimal() + scale_colour_manual(values = c("mediumaquamarine", "maroon4")) +
	theme(plot.title = element_text(family="Times", face="bold", size=20)) +
	ylab("Surface Area (mm2)") +
	labs(title = "Parietal Surface Area (Female)", subtitle = mainstring)

# Temporal
stats_temporal_cort_F <- lmer(Temporal_CT ~ Age*psTerminal + (Age | bblid), df[df$sex == 2,])
coeffdf <- summary(stats_temporal_cort_F)$coefficients
estimates <- coeffdf[,"Estimate"]
stderrors <- coeffdf[,"Std. Error"]
tvalues <- coeffdf[,"t value"]
if (abs(tvalues[["Age"]]) > 2) { age_est <- paste0(round(estimates[["Age"]], digits=3), "*") } else { age_est <- round(estimates[["Age"]], digits=3) }
if (abs(tvalues[["psTerminalTD"]]) > 2) { diag_est <- paste0(round(estimates[["psTerminalTD"]], digits=3), "*") } else { diag_est <- round(estimates[["psTerminalTD"]], digits=3) }
if (abs(tvalues[["Age:psTerminalTD"]]) > 2) { inter_est <- paste0(round(estimates[["Age:psTerminalTD"]], digits=3), "*") } else { inter_est <- round(estimates[["Age:psTerminalTD"]], digits=3) }
mainstring <- paste0(femalestring, "\n", "Coefficients: Age=", age_est, ", Diagnosis=", diag_est, ", Age:Diagnosis=", inter_est)
temporal_cort_F <- ggplot(df[df$sex == 2,], aes(x=Age, y=Temporal_CT)) + 
	geom_line(aes(group=bblid, col=psTerminal), alpha=.2) + 
	geom_smooth(se = TRUE, method = "glm", aes(group=psTerminal, col=psTerminal)) + 
	theme_minimal() + scale_colour_manual(values = c("mediumaquamarine", "maroon4")) +
	theme(plot.title = element_text(family="Times", face="bold", size=20)) +
	ylab("Cortical Thickness (mm)") +
	labs(title = "Temporal Cortical Thickness (Female)", subtitle = mainstring)

stats_temporal_vol_F <- lmer(Temporal_Vol ~ Age*psTerminal + (Age | bblid), df[df$sex == 2,])
coeffdf <- summary(stats_temporal_vol_F)$coefficients
estimates <- coeffdf[,"Estimate"]
stderrors <- coeffdf[,"Std. Error"]
tvalues <- coeffdf[,"t value"]
if (abs(tvalues[["Age"]]) > 2) { age_est <- paste0(round(estimates[["Age"]], digits=3), "*") } else { age_est <- round(estimates[["Age"]], digits=3) }
if (abs(tvalues[["psTerminalTD"]]) > 2) { diag_est <- paste0(round(estimates[["psTerminalTD"]], digits=3), "*") } else { diag_est <- round(estimates[["psTerminalTD"]], digits=3) }
if (abs(tvalues[["Age:psTerminalTD"]]) > 2) { inter_est <- paste0(round(estimates[["Age:psTerminalTD"]], digits=3), "*") } else { inter_est <- round(estimates[["Age:psTerminalTD"]], digits=3) }
mainstring <- paste0(femalestring, "\n", "Coefficients: Age=", age_est, ", Diagnosis=", diag_est, ", Age:Diagnosis=", inter_est)
temporal_vol_F <- ggplot(df[df$sex == 2,], aes(x=Age, y=Temporal_Vol)) + 
	geom_line(aes(group=bblid, col=psTerminal), alpha=.2) + 
	geom_smooth(se = TRUE, method = "glm", aes(group=psTerminal, col=psTerminal)) + 
	theme_minimal() + scale_colour_manual(values = c("mediumaquamarine", "maroon4")) +
	theme(plot.title = element_text(family="Times", face="bold", size=20)) +
	ylab("Gray Matter Volume (mm3)") +
	labs(title = "Temporal Gray Matter Volume (Female)", subtitle = mainstring)

stats_temporal_gmd_F <- lmer(Temporal_GMD ~ Age*psTerminal + (Age | bblid), df[df$sex == 2,])
coeffdf <- summary(stats_temporal_gmd_F)$coefficients
estimates <- coeffdf[,"Estimate"]
stderrors <- coeffdf[,"Std. Error"]
tvalues <- coeffdf[,"t value"]
if (abs(tvalues[["Age"]]) > 2) { age_est <- paste0(round(estimates[["Age"]], digits=3), "*") } else { age_est <- round(estimates[["Age"]], digits=3) }
if (abs(tvalues[["psTerminalTD"]]) > 2) { diag_est <- paste0(round(estimates[["psTerminalTD"]], digits=3), "*") } else { diag_est <- round(estimates[["psTerminalTD"]], digits=3) }
if (abs(tvalues[["Age:psTerminalTD"]]) > 2) { inter_est <- paste0(round(estimates[["Age:psTerminalTD"]], digits=3), "*") } else { inter_est <- round(estimates[["Age:psTerminalTD"]], digits=3) }
mainstring <- paste0(femalestring, "\n", "Coefficients: Age=", age_est, ", Diagnosis=", diag_est, ", Age:Diagnosis=", inter_est)
temporal_gmd_F <- ggplot(df[df$sex == 2,], aes(x=Age, y=Temporal_GMD)) + 
	geom_line(aes(group=bblid, col=psTerminal), alpha=.2) + 
	geom_smooth(se = TRUE, method = "glm", aes(group=psTerminal, col=psTerminal)) + 
	theme_minimal() + scale_colour_manual(values = c("mediumaquamarine", "maroon4")) +
	theme(plot.title = element_text(family="Times", face="bold", size=20)) +
	ylab("Gray Matter Density") +
	labs(title = "Temporal Gray Matter Density (Female)", subtitle = mainstring)

stats_temporal_sa_F <- lmer(Temporal_SA ~ Age*psTerminal + (Age | bblid), df[df$sex == 2,])
coeffdf <- summary(stats_temporal_sa_F)$coefficients
estimates <- coeffdf[,"Estimate"]
stderrors <- coeffdf[,"Std. Error"]
tvalues <- coeffdf[,"t value"]
if (abs(tvalues[["Age"]]) > 2) { age_est <- paste0(round(estimates[["Age"]], digits=3), "*") } else { age_est <- round(estimates[["Age"]], digits=3) }
if (abs(tvalues[["psTerminalTD"]]) > 2) { diag_est <- paste0(round(estimates[["psTerminalTD"]], digits=3), "*") } else { diag_est <- round(estimates[["psTerminalTD"]], digits=3) }
if (abs(tvalues[["Age:psTerminalTD"]]) > 2) { inter_est <- paste0(round(estimates[["Age:psTerminalTD"]], digits=3), "*") } else { inter_est <- round(estimates[["Age:psTerminalTD"]], digits=3) }
mainstring <- paste0(femalestring, "\n", "Coefficients: Age=", age_est, ", Diagnosis=", diag_est, ", Age:Diagnosis=", inter_est)
temporal_sa_F <- ggplot(df[df$sex == 2,], aes(x=Age, y=Temporal_SA)) + 
	geom_line(aes(group=bblid, col=psTerminal), alpha=.2) + 
	geom_smooth(se = TRUE, method = "glm", aes(group=psTerminal, col=psTerminal)) + 
	theme_minimal() + scale_colour_manual(values = c("mediumaquamarine", "maroon4")) +
	theme(plot.title = element_text(family="Times", face="bold", size=20)) +
	ylab("Surface Area (mm2)") +
	labs(title = "Temporal Surface Area (Female)", subtitle = mainstring)

# Occipital
stats_occipital_cort_F <- lmer(Occipital_CT ~ Age*psTerminal + (Age | bblid), df[df$sex == 2,])
coeffdf <- summary(stats_occipital_cort_F)$coefficients
estimates <- coeffdf[,"Estimate"]
stderrors <- coeffdf[,"Std. Error"]
tvalues <- coeffdf[,"t value"]
if (abs(tvalues[["Age"]]) > 2) { age_est <- paste0(round(estimates[["Age"]], digits=3), "*") } else { age_est <- round(estimates[["Age"]], digits=3) }
if (abs(tvalues[["psTerminalTD"]]) > 2) { diag_est <- paste0(round(estimates[["psTerminalTD"]], digits=3), "*") } else { diag_est <- round(estimates[["psTerminalTD"]], digits=3) }
if (abs(tvalues[["Age:psTerminalTD"]]) > 2) { inter_est <- paste0(round(estimates[["Age:psTerminalTD"]], digits=3), "*") } else { inter_est <- round(estimates[["Age:psTerminalTD"]], digits=3) }
mainstring <- paste0(femalestring, "\n", "Coefficients: Age=", age_est, ", Diagnosis=", diag_est, ", Age:Diagnosis=", inter_est)
occipital_cort_F <- ggplot(df[df$sex == 2,], aes(x=Age, y=Occipital_CT)) + 
	geom_line(aes(group=bblid, col=psTerminal), alpha=.2) + 
	geom_smooth(se = TRUE, method = "glm", aes(group=psTerminal, col=psTerminal)) + 
	theme_minimal() + scale_colour_manual(values = c("mediumaquamarine", "maroon4")) +
	theme(plot.title = element_text(family="Times", face="bold", size=20)) +
	ylab("Cortical Thickness (mm)") +
	labs(title = "Occipital Cortical Thickness (Female)", subtitle = mainstring)

stats_occipital_vol_F <- lmer(Occipital_Vol ~ Age*psTerminal + (Age | bblid), df[df$sex == 2,])
coeffdf <- summary(stats_occipital_vol_F)$coefficients
estimates <- coeffdf[,"Estimate"]
stderrors <- coeffdf[,"Std. Error"]
tvalues <- coeffdf[,"t value"]
if (abs(tvalues[["Age"]]) > 2) { age_est <- paste0(round(estimates[["Age"]], digits=3), "*") } else { age_est <- round(estimates[["Age"]], digits=3) }
if (abs(tvalues[["psTerminalTD"]]) > 2) { diag_est <- paste0(round(estimates[["psTerminalTD"]], digits=3), "*") } else { diag_est <- round(estimates[["psTerminalTD"]], digits=3) }
if (abs(tvalues[["Age:psTerminalTD"]]) > 2) { inter_est <- paste0(round(estimates[["Age:psTerminalTD"]], digits=3), "*") } else { inter_est <- round(estimates[["Age:psTerminalTD"]], digits=3) }
mainstring <- paste0(femalestring, "\n", "Coefficients: Age=", age_est, ", Diagnosis=", diag_est, ", Age:Diagnosis=", inter_est)
occipital_vol_F <- ggplot(df[df$sex == 2,], aes(x=Age, y=Occipital_Vol)) + 
	geom_line(aes(group=bblid, col=psTerminal), alpha=.2) + 
	geom_smooth(se = TRUE, method = "glm", aes(group=psTerminal, col=psTerminal)) + 
	theme_minimal() + scale_colour_manual(values = c("mediumaquamarine", "maroon4")) +
	theme(plot.title = element_text(family="Times", face="bold", size=20)) +
	ylab("Gray Matter Volume (mm3)") +
	labs(title = "Occipital Gray Matter Volume (Female)", subtitle = mainstring)

stats_occipital_gmd_F <- lmer(Occipital_GMD ~ Age*psTerminal + (Age | bblid), df[df$sex == 2,])
coeffdf <- summary(stats_occipital_gmd_F)$coefficients
estimates <- coeffdf[,"Estimate"]
stderrors <- coeffdf[,"Std. Error"]
tvalues <- coeffdf[,"t value"]
if (abs(tvalues[["Age"]]) > 2) { age_est <- paste0(round(estimates[["Age"]], digits=3), "*") } else { age_est <- round(estimates[["Age"]], digits=3) }
if (abs(tvalues[["psTerminalTD"]]) > 2) { diag_est <- paste0(round(estimates[["psTerminalTD"]], digits=3), "*") } else { diag_est <- round(estimates[["psTerminalTD"]], digits=3) }
if (abs(tvalues[["Age:psTerminalTD"]]) > 2) { inter_est <- paste0(round(estimates[["Age:psTerminalTD"]], digits=3), "*") } else { inter_est <- round(estimates[["Age:psTerminalTD"]], digits=3) }
mainstring <- paste0(femalestring, "\n", "Coefficients: Age=", age_est, ", Diagnosis=", diag_est, ", Age:Diagnosis=", inter_est)
occipital_gmd_F <- ggplot(df[df$sex == 2,], aes(x=Age, y=Occipital_GMD)) + 
	geom_line(aes(group=bblid, col=psTerminal), alpha=.2) + 
	geom_smooth(se = TRUE, method = "glm", aes(group=psTerminal, col=psTerminal)) + 
	theme_minimal() + scale_colour_manual(values = c("mediumaquamarine", "maroon4")) +
	theme(plot.title = element_text(family="Times", face="bold", size=20)) +
	ylab("Gray Matter Density") +
	labs(title = "Occipital Gray Matter Density (Female)", subtitle = mainstring)

stats_occipital_sa_F <- lmer(Occipital_SA ~ Age*psTerminal + (Age | bblid), df[df$sex == 2,])
coeffdf <- summary(stats_occipital_sa_F)$coefficients
estimates <- coeffdf[,"Estimate"]
stderrors <- coeffdf[,"Std. Error"]
tvalues <- coeffdf[,"t value"]
if (abs(tvalues[["Age"]]) > 2) { age_est <- paste0(round(estimates[["Age"]], digits=3), "*") } else { age_est <- round(estimates[["Age"]], digits=3) }
if (abs(tvalues[["psTerminalTD"]]) > 2) { diag_est <- paste0(round(estimates[["psTerminalTD"]], digits=3), "*") } else { diag_est <- round(estimates[["psTerminalTD"]], digits=3) }
if (abs(tvalues[["Age:psTerminalTD"]]) > 2) { inter_est <- paste0(round(estimates[["Age:psTerminalTD"]], digits=3), "*") } else { inter_est <- round(estimates[["Age:psTerminalTD"]], digits=3) }
mainstring <- paste0(femalestring, "\n", "Coefficients: Age=", age_est, ", Diagnosis=", diag_est, ", Age:Diagnosis=", inter_est)
occipital_sa_F <- ggplot(df[df$sex == 2,], aes(x=Age, y=Occipital_SA)) + 
	geom_line(aes(group=bblid, col=psTerminal), alpha=.2) + 
	geom_smooth(se = TRUE, method = "glm", aes(group=psTerminal, col=psTerminal)) + 
	theme_minimal() + scale_colour_manual(values = c("mediumaquamarine", "maroon4")) +
	theme(plot.title = element_text(family="Times", face="bold", size=20)) +
	ylab("Surface Area (mm2)") +
	labs(title = "Occipital Surface Area (Female)", subtitle = mainstring)

# BasGang
stats_basGang_vol_F <- lmer(BasGang_Vol ~ Age*psTerminal + (Age | bblid), df[df$sex == 2,])
coeffdf <- summary(stats_basGang_vol_F)$coefficients
estimates <- coeffdf[,"Estimate"]
stderrors <- coeffdf[,"Std. Error"]
tvalues <- coeffdf[,"t value"]
if (abs(tvalues[["Age"]]) > 2) { age_est <- paste0(round(estimates[["Age"]], digits=3), "*") } else { age_est <- round(estimates[["Age"]], digits=3) }
if (abs(tvalues[["psTerminalTD"]]) > 2) { diag_est <- paste0(round(estimates[["psTerminalTD"]], digits=3), "*") } else { diag_est <- round(estimates[["psTerminalTD"]], digits=3) }
if (abs(tvalues[["Age:psTerminalTD"]]) > 2) { inter_est <- paste0(round(estimates[["Age:psTerminalTD"]], digits=3), "*") } else { inter_est <- round(estimates[["Age:psTerminalTD"]], digits=3) }
mainstring <- paste0(femalestring, "\n", "Coefficients: Age=", age_est, ", Diagnosis=", diag_est, ", Age:Diagnosis=", inter_est)
basGang_vol_F <- ggplot(df[df$sex == 2,], aes(x=Age, y=BasGang_Vol)) + 
	geom_line(aes(group=bblid, col=psTerminal), alpha=.2) + 
	geom_smooth(se = TRUE, method = "glm", aes(group=psTerminal, col=psTerminal)) + 
	theme_minimal() + scale_colour_manual(values = c("mediumaquamarine", "maroon4")) +
	theme(plot.title = element_text(family="Times", face="bold", size=20)) +
	ylab("Gray Matter Volume (mm3)") +
	labs(title = "Basal Gangalia Gray Matter Volume (Female)", subtitle = mainstring)

stats_basGang_gmd_F <- lmer(BasGang_GMD ~ Age*psTerminal + (Age | bblid), df[df$sex == 2,])
coeffdf <- summary(stats_basGang_gmd_F)$coefficients
estimates <- coeffdf[,"Estimate"]
stderrors <- coeffdf[,"Std. Error"]
tvalues <- coeffdf[,"t value"]
if (abs(tvalues[["Age"]]) > 2) { age_est <- paste0(round(estimates[["Age"]], digits=3), "*") } else { age_est <- round(estimates[["Age"]], digits=3) }
if (abs(tvalues[["psTerminalTD"]]) > 2) { diag_est <- paste0(round(estimates[["psTerminalTD"]], digits=3), "*") } else { diag_est <- round(estimates[["psTerminalTD"]], digits=3) }
if (abs(tvalues[["Age:psTerminalTD"]]) > 2) { inter_est <- paste0(round(estimates[["Age:psTerminalTD"]], digits=3), "*") } else { inter_est <- round(estimates[["Age:psTerminalTD"]], digits=3) }
mainstring <- paste0(femalestring, "\n", "Coefficients: Age=", age_est, ", Diagnosis=", diag_est, ", Age:Diagnosis=", inter_est)
basGang_gmd_F <- ggplot(df[df$sex == 2,], aes(x=Age, y=BasGang_GMD)) + 
	geom_line(aes(group=bblid, col=psTerminal), alpha=.2) + 
	geom_smooth(se = TRUE, method = "glm", aes(group=psTerminal, col=psTerminal)) + 
	theme_minimal() + scale_colour_manual(values = c("mediumaquamarine", "maroon4")) +
	theme(plot.title = element_text(family="Times", face="bold", size=20)) +
	ylab("Gray Matter Density") +
	labs(title = "Basal Gangalia Gray Matter Density (Female)", subtitle = mainstring)

# Limbic
stats_limbic_vol_F <- lmer(Limbic_Vol ~ Age*psTerminal + (Age | bblid), df[df$sex == 2,])
coeffdf <- summary(stats_limbic_vol_F)$coefficients
estimates <- coeffdf[,"Estimate"]
stderrors <- coeffdf[,"Std. Error"]
tvalues <- coeffdf[,"t value"]
if (abs(tvalues[["Age"]]) > 2) { age_est <- paste0(round(estimates[["Age"]], digits=3), "*") } else { age_est <- round(estimates[["Age"]], digits=3) }
if (abs(tvalues[["psTerminalTD"]]) > 2) { diag_est <- paste0(round(estimates[["psTerminalTD"]], digits=3), "*") } else { diag_est <- round(estimates[["psTerminalTD"]], digits=3) }
if (abs(tvalues[["Age:psTerminalTD"]]) > 2) { inter_est <- paste0(round(estimates[["Age:psTerminalTD"]], digits=3), "*") } else { inter_est <- round(estimates[["Age:psTerminalTD"]], digits=3) }
mainstring <- paste0(femalestring, "\n", "Coefficients: Age=", age_est, ", Diagnosis=", diag_est, ", Age:Diagnosis=", inter_est)
limbic_vol_F <- ggplot(df[df$sex == 2,], aes(x=Age, y=Limbic_Vol)) + 
	geom_line(aes(group=bblid, col=psTerminal), alpha=.2) + 
	geom_smooth(se = TRUE, method = "glm", aes(group=psTerminal, col=psTerminal)) + 
	theme_minimal() + scale_colour_manual(values = c("mediumaquamarine", "maroon4")) +
	theme(plot.title = element_text(family="Times", face="bold", size=20)) +
	ylab("Gray Matter Volume (mm3)") +
	labs(title = "Limbic Gray Matter Volume (Female)", subtitle = mainstring)

stats_limbic_gmd_F <- lmer(Limbic_GMD ~ Age*psTerminal + (Age | bblid), df[df$sex == 2,])
coeffdf <- summary(stats_limbic_gmd_F)$coefficients
estimates <- coeffdf[,"Estimate"]
stderrors <- coeffdf[,"Std. Error"]
tvalues <- coeffdf[,"t value"]
if (abs(tvalues[["Age"]]) > 2) { age_est <- paste0(round(estimates[["Age"]], digits=3), "*") } else { age_est <- round(estimates[["Age"]], digits=3) }
if (abs(tvalues[["psTerminalTD"]]) > 2) { diag_est <- paste0(round(estimates[["psTerminalTD"]], digits=3), "*") } else { diag_est <- round(estimates[["psTerminalTD"]], digits=3) }
if (abs(tvalues[["Age:psTerminalTD"]]) > 2) { inter_est <- paste0(round(estimates[["Age:psTerminalTD"]], digits=3), "*") } else { inter_est <- round(estimates[["Age:psTerminalTD"]], digits=3) }
mainstring <- paste0(femalestring, "\n", "Coefficients: Age=", age_est, ", Diagnosis=", diag_est, ", Age:Diagnosis=", inter_est)
limbic_gmd_F <- ggplot(df[df$sex == 2,], aes(x=Age, y=Limbic_GMD)) + 
	geom_line(aes(group=bblid, col=psTerminal), alpha=.2) + 
	geom_smooth(se = TRUE, method = "glm", aes(group=psTerminal, col=psTerminal)) + 
	theme_minimal() + scale_colour_manual(values = c("mediumaquamarine", "maroon4")) +
	theme(plot.title = element_text(family="Times", face="bold", size=20)) +
	ylab("Gray Matter Density") +
	labs(title = "Limbic Gray Matter Density (Female)", subtitle = mainstring)

#----------------------------------- Male -----------------------------------#
# Frontal Plots
stats_frontorb_cort_M <- lmer(FrontOrb_CT ~ Age*psTerminal + (Age | bblid), df[df$sex == 2,])
coeffdf <- summary(stats_frontorb_cort_M)$coefficients
estimates <- coeffdf[,"Estimate"]
stderrors <- coeffdf[,"Std. Error"]
tvalues <- coeffdf[,"t value"]
if (abs(tvalues[["Age"]]) > 2) { age_est <- paste0(round(estimates[["Age"]], digits=3), "*") } else { age_est <- round(estimates[["Age"]], digits=3) }
if (abs(tvalues[["psTerminalTD"]]) > 2) { diag_est <- paste0(round(estimates[["psTerminalTD"]], digits=3), "*") } else { diag_est <- round(estimates[["psTerminalTD"]], digits=3) }
if (abs(tvalues[["Age:psTerminalTD"]]) > 2) { inter_est <- paste0(round(estimates[["Age:psTerminalTD"]], digits=3), "*") } else { inter_est <- round(estimates[["Age:psTerminalTD"]], digits=3) }
mainstring <- paste0(malestring, "\n", "Coefficients: Age=", age_est, ", Diagnosis=", diag_est, ", Age:Diagnosis=", inter_est)
frontorb_cort_M <- ggplot(df[df$sex == 2,], aes(x=Age, y=FrontOrb_CT)) + 
	geom_line(aes(group=bblid, col=psTerminal), alpha=.2) + 
	geom_smooth(se = TRUE, method = "glm", aes(group=psTerminal, col=psTerminal)) + 
	theme_minimal() + scale_colour_manual(values = c("mediumaquamarine", "maroon4")) +
	theme(plot.title = element_text(family="Times", face="bold", size=20)) +
	ylab("Cortical Thickness (mm)") +
	labs(title = "Frontal Orbital Cortical Thickness (Male)", subtitle = mainstring)

stats_frontorb_vol_M <- lmer(FrontOrb_Vol ~ Age*psTerminal + (Age | bblid), df[df$sex == 2,])
coeffdf <- summary(stats_frontorb_vol_M)$coefficients
estimates <- coeffdf[,"Estimate"]
stderrors <- coeffdf[,"Std. Error"]
tvalues <- coeffdf[,"t value"]
if (abs(tvalues[["Age"]]) > 2) { age_est <- paste0(round(estimates[["Age"]], digits=3), "*") } else { age_est <- round(estimates[["Age"]], digits=3) }
if (abs(tvalues[["psTerminalTD"]]) > 2) { diag_est <- paste0(round(estimates[["psTerminalTD"]], digits=3), "*") } else { diag_est <- round(estimates[["psTerminalTD"]], digits=3) }
if (abs(tvalues[["Age:psTerminalTD"]]) > 2) { inter_est <- paste0(round(estimates[["Age:psTerminalTD"]], digits=3), "*") } else { inter_est <- round(estimates[["Age:psTerminalTD"]], digits=3) }
mainstring <- paste0(malestring, "\n", "Coefficients: Age=", age_est, ", Diagnosis=", diag_est, ", Age:Diagnosis=", inter_est)
frontorb_vol_M <- ggplot(df[df$sex == 2,], aes(x=Age, y=FrontOrb_Vol)) + 
	geom_line(aes(group=bblid, col=psTerminal), alpha=.2) + 
	geom_smooth(se = TRUE, method = "glm", aes(group=psTerminal, col=psTerminal)) + 
	theme_minimal() + scale_colour_manual(values = c("mediumaquamarine", "maroon4")) +
	theme(plot.title = element_text(family="Times", face="bold", size=20)) +
	ylab("Gray Matter Volume (mm3)") +
	labs(title = "Frontal Orbital Gray Matter Volume (Male)", subtitle = mainstring)

stats_frontorb_gmd_M <- lmer(FrontOrb_GMD ~ Age*psTerminal + (Age | bblid), df[df$sex == 2,])
coeffdf <- summary(stats_frontorb_gmd_M)$coefficients
estimates <- coeffdf[,"Estimate"]
stderrors <- coeffdf[,"Std. Error"]
tvalues <- coeffdf[,"t value"]
if (abs(tvalues[["Age"]]) > 2) { age_est <- paste0(round(estimates[["Age"]], digits=3), "*") } else { age_est <- round(estimates[["Age"]], digits=3) }
if (abs(tvalues[["psTerminalTD"]]) > 2) { diag_est <- paste0(round(estimates[["psTerminalTD"]], digits=3), "*") } else { diag_est <- round(estimates[["psTerminalTD"]], digits=3) }
if (abs(tvalues[["Age:psTerminalTD"]]) > 2) { inter_est <- paste0(round(estimates[["Age:psTerminalTD"]], digits=3), "*") } else { inter_est <- round(estimates[["Age:psTerminalTD"]], digits=3) }
mainstring <- paste0(malestring, "\n", "Coefficients: Age=", age_est, ", Diagnosis=", diag_est, ", Age:Diagnosis=", inter_est)
frontorb_gmd_M <- ggplot(df[df$sex == 2,], aes(x=Age, y=FrontOrb_GMD)) + 
	geom_line(aes(group=bblid, col=psTerminal), alpha=.2) + 
	geom_smooth(se = TRUE, method = "glm", aes(group=psTerminal, col=psTerminal)) + 
	theme_minimal() + scale_colour_manual(values = c("mediumaquamarine", "maroon4")) +
	theme(plot.title = element_text(family="Times", face="bold", size=20)) +
	ylab("Gray Matter Density") +
	labs(title = "Frontal Orbital Gray Matter Density (Male)", subtitle = mainstring)

stats_frontdors_cort_M <- lmer(FrontDors_CT ~ Age*psTerminal + (Age | bblid), df[df$sex == 2,])
coeffdf <- summary(stats_frontdors_cort_M)$coefficients
estimates <- coeffdf[,"Estimate"]
stderrors <- coeffdf[,"Std. Error"]
tvalues <- coeffdf[,"t value"]
if (abs(tvalues[["Age"]]) > 2) { age_est <- paste0(round(estimates[["Age"]], digits=3), "*") } else { age_est <- round(estimates[["Age"]], digits=3) }
if (abs(tvalues[["psTerminalTD"]]) > 2) { diag_est <- paste0(round(estimates[["psTerminalTD"]], digits=3), "*") } else { diag_est <- round(estimates[["psTerminalTD"]], digits=3) }
if (abs(tvalues[["Age:psTerminalTD"]]) > 2) { inter_est <- paste0(round(estimates[["Age:psTerminalTD"]], digits=3), "*") } else { inter_est <- round(estimates[["Age:psTerminalTD"]], digits=3) }
mainstring <- paste0(malestring, "\n", "Coefficients: Age=", age_est, ", Diagnosis=", diag_est, ", Age:Diagnosis=", inter_est)
frontdors_cort_M <- ggplot(df[df$sex == 2,], aes(x=Age, y=FrontDors_CT)) + 
	geom_line(aes(group=bblid, col=psTerminal), alpha=.2) + 
	geom_smooth(se = TRUE, method = "glm", aes(group=psTerminal, col=psTerminal)) + 
	theme_minimal() + scale_colour_manual(values = c("mediumaquamarine", "maroon4")) +
	theme(plot.title = element_text(family="Times", face="bold", size=20)) +
	ylab("Cortical Thickness (mm)") +
	labs(title = "Frontal Dorsal Cortical Thickness (Male)", subtitle = mainstring)

stats_frontdors_vol_M <- lmer(FrontDors_Vol ~ Age*psTerminal + (Age | bblid), df[df$sex == 2,])
coeffdf <- summary(stats_frontdors_vol_M)$coefficients
estimates <- coeffdf[,"Estimate"]
stderrors <- coeffdf[,"Std. Error"]
tvalues <- coeffdf[,"t value"]
if (abs(tvalues[["Age"]]) > 2) { age_est <- paste0(round(estimates[["Age"]], digits=3), "*") } else { age_est <- round(estimates[["Age"]], digits=3) }
if (abs(tvalues[["psTerminalTD"]]) > 2) { diag_est <- paste0(round(estimates[["psTerminalTD"]], digits=3), "*") } else { diag_est <- round(estimates[["psTerminalTD"]], digits=3) }
if (abs(tvalues[["Age:psTerminalTD"]]) > 2) { inter_est <- paste0(round(estimates[["Age:psTerminalTD"]], digits=3), "*") } else { inter_est <- round(estimates[["Age:psTerminalTD"]], digits=3) }
mainstring <- paste0(malestring, "\n", "Coefficients: Age=", age_est, ", Diagnosis=", diag_est, ", Age:Diagnosis=", inter_est)
frontdors_vol_M <- ggplot(df[df$sex == 2,], aes(x=Age, y=FrontDors_Vol)) + 
	geom_line(aes(group=bblid, col=psTerminal), alpha=.2) + 
	geom_smooth(se = TRUE, method = "glm", aes(group=psTerminal, col=psTerminal)) + 
	theme_minimal() + scale_colour_manual(values = c("mediumaquamarine", "maroon4")) +
	theme(plot.title = element_text(family="Times", face="bold", size=20)) +
	ylab("Gray Matter Volume (mm3)") +
	labs(title = "Frontal Dorsal Gray Matter Volume (Male)", subtitle = mainstring)

stats_frontdors_gmd_M <- lmer(FrontDors_GMD ~ Age*psTerminal + (Age | bblid), df[df$sex == 2,])
coeffdf <- summary(stats_frontdors_gmd_M)$coefficients
estimates <- coeffdf[,"Estimate"]
stderrors <- coeffdf[,"Std. Error"]
tvalues <- coeffdf[,"t value"]
if (abs(tvalues[["Age"]]) > 2) { age_est <- paste0(round(estimates[["Age"]], digits=3), "*") } else { age_est <- round(estimates[["Age"]], digits=3) }
if (abs(tvalues[["psTerminalTD"]]) > 2) { diag_est <- paste0(round(estimates[["psTerminalTD"]], digits=3), "*") } else { diag_est <- round(estimates[["psTerminalTD"]], digits=3) }
if (abs(tvalues[["Age:psTerminalTD"]]) > 2) { inter_est <- paste0(round(estimates[["Age:psTerminalTD"]], digits=3), "*") } else { inter_est <- round(estimates[["Age:psTerminalTD"]], digits=3) }
mainstring <- paste0(malestring, "\n", "Coefficients: Age=", age_est, ", Diagnosis=", diag_est, ", Age:Diagnosis=", inter_est)
frontdors_gmd_M <- ggplot(df[df$sex == 2,], aes(x=Age, y=FrontDors_GMD)) + 
	geom_line(aes(group=bblid, col=psTerminal), alpha=.2) + 
	geom_smooth(se = TRUE, method = "glm", aes(group=psTerminal, col=psTerminal)) + 
	theme_minimal() + scale_colour_manual(values = c("mediumaquamarine", "maroon4")) +
	theme(plot.title = element_text(family="Times", face="bold", size=20)) +
	ylab("Gray Matter Density") +
	labs(title = "Frontal Dorsal Gray Matter Density (Male)", subtitle = mainstring)

stats_front_sa_M <- lmer(Frontal_SA ~ Age*psTerminal + (Age | bblid), df[df$sex == 2,])
coeffdf <- summary(stats_front_sa_M)$coefficients
estimates <- coeffdf[,"Estimate"]
stderrors <- coeffdf[,"Std. Error"]
tvalues <- coeffdf[,"t value"]
if (abs(tvalues[["Age"]]) > 2) { age_est <- paste0(round(estimates[["Age"]], digits=3), "*") } else { age_est <- round(estimates[["Age"]], digits=3) }
if (abs(tvalues[["psTerminalTD"]]) > 2) { diag_est <- paste0(round(estimates[["psTerminalTD"]], digits=3), "*") } else { diag_est <- round(estimates[["psTerminalTD"]], digits=3) }
if (abs(tvalues[["Age:psTerminalTD"]]) > 2) { inter_est <- paste0(round(estimates[["Age:psTerminalTD"]], digits=3), "*") } else { inter_est <- round(estimates[["Age:psTerminalTD"]], digits=3) }
mainstring <- paste0(malestring, "\n", "Coefficients: Age=", age_est, ", Diagnosis=", diag_est, ", Age:Diagnosis=", inter_est)
front_sa_M <- ggplot(df[df$sex == 2,], aes(x=Age, y=Frontal_SA)) + 
	geom_line(aes(group=bblid, col=psTerminal), alpha=.2) + 
	geom_smooth(se = TRUE, method = "glm", aes(group=psTerminal, col=psTerminal)) + 
	theme_minimal() + scale_colour_manual(values = c("mediumaquamarine", "maroon4")) +
	theme(plot.title = element_text(family="Times", face="bold", size=20)) +
	ylab("Surface Area (mm2)") +
	labs(title = "Frontal Surface Area (Male)", subtitle = mainstring)

# Parietal
stats_parietal_cort_M <- lmer(Parietal_CT ~ Age*psTerminal + (Age | bblid), df[df$sex == 2,])
coeffdf <- summary(stats_parietal_cort_M)$coefficients
estimates <- coeffdf[,"Estimate"]
stderrors <- coeffdf[,"Std. Error"]
tvalues <- coeffdf[,"t value"]
if (abs(tvalues[["Age"]]) > 2) { age_est <- paste0(round(estimates[["Age"]], digits=3), "*") } else { age_est <- round(estimates[["Age"]], digits=3) }
if (abs(tvalues[["psTerminalTD"]]) > 2) { diag_est <- paste0(round(estimates[["psTerminalTD"]], digits=3), "*") } else { diag_est <- round(estimates[["psTerminalTD"]], digits=3) }
if (abs(tvalues[["Age:psTerminalTD"]]) > 2) { inter_est <- paste0(round(estimates[["Age:psTerminalTD"]], digits=3), "*") } else { inter_est <- round(estimates[["Age:psTerminalTD"]], digits=3) }
mainstring <- paste0(malestring, "\n", "Coefficients: Age=", age_est, ", Diagnosis=", diag_est, ", Age:Diagnosis=", inter_est)
parietal_cort_M <- ggplot(df[df$sex == 2,], aes(x=Age, y=Parietal_CT)) + 
	geom_line(aes(group=bblid, col=psTerminal), alpha=.2) + 
	geom_smooth(se = TRUE, method = "glm",  aes(group=psTerminal, col=psTerminal)) + 
	theme_minimal() + scale_colour_manual(values = c("mediumaquamarine", "maroon4")) +
	theme(plot.title = element_text(family="Times", face="bold", size=20)) +
	ylab("Cortical Thickness (mm)") +
	labs(title = "Parietal Cortical Thickness (Male)", subtitle = mainstring)

stats_parietal_vol_M <- lmer(Parietal_Vol ~ Age*psTerminal + (Age | bblid), df[df$sex == 2,])
coeffdf <- summary(stats_parietal_vol_M)$coefficients
estimates <- coeffdf[,"Estimate"]
stderrors <- coeffdf[,"Std. Error"]
tvalues <- coeffdf[,"t value"]
if (abs(tvalues[["Age"]]) > 2) { age_est <- paste0(round(estimates[["Age"]], digits=3), "*") } else { age_est <- round(estimates[["Age"]], digits=3) }
if (abs(tvalues[["psTerminalTD"]]) > 2) { diag_est <- paste0(round(estimates[["psTerminalTD"]], digits=3), "*") } else { diag_est <- round(estimates[["psTerminalTD"]], digits=3) }
if (abs(tvalues[["Age:psTerminalTD"]]) > 2) { inter_est <- paste0(round(estimates[["Age:psTerminalTD"]], digits=3), "*") } else { inter_est <- round(estimates[["Age:psTerminalTD"]], digits=3) }
mainstring <- paste0(malestring, "\n", "Coefficients: Age=", age_est, ", Diagnosis=", diag_est, ", Age:Diagnosis=", inter_est)
parietal_vol_M <- ggplot(df[df$sex == 2,], aes(x=Age, y=Parietal_Vol)) + 
	geom_line(aes(group=bblid, col=psTerminal), alpha=.2) + 
	geom_smooth(se = TRUE, method = "glm", aes(group=psTerminal, col=psTerminal)) + 
	theme_minimal() + scale_colour_manual(values = c("mediumaquamarine", "maroon4")) +
	theme(plot.title = element_text(family="Times", face="bold", size=20)) +
	ylab("Gray Matter Volume (mm3)") +
	labs(title = "Parietal Gray Matter Volume (Male)", subtitle = mainstring)

stats_parietal_gmd_M <- lmer(Parietal_GMD ~ Age*psTerminal + (Age | bblid), df[df$sex == 2,])
coeffdf <- summary(stats_parietal_gmd_M)$coefficients
estimates <- coeffdf[,"Estimate"]
stderrors <- coeffdf[,"Std. Error"]
tvalues <- coeffdf[,"t value"]
if (abs(tvalues[["Age"]]) > 2) { age_est <- paste0(round(estimates[["Age"]], digits=3), "*") } else { age_est <- round(estimates[["Age"]], digits=3) }
if (abs(tvalues[["psTerminalTD"]]) > 2) { diag_est <- paste0(round(estimates[["psTerminalTD"]], digits=3), "*") } else { diag_est <- round(estimates[["psTerminalTD"]], digits=3) }
if (abs(tvalues[["Age:psTerminalTD"]]) > 2) { inter_est <- paste0(round(estimates[["Age:psTerminalTD"]], digits=3), "*") } else { inter_est <- round(estimates[["Age:psTerminalTD"]], digits=3) }
mainstring <- paste0(malestring, "\n", "Coefficients: Age=", age_est, ", Diagnosis=", diag_est, ", Age:Diagnosis=", inter_est)
parietal_gmd_M <- ggplot(df[df$sex == 2,], aes(x=Age, y=Parietal_GMD)) + 
	geom_line(aes(group=bblid, col=psTerminal), alpha=.2) + 
	geom_smooth(se = TRUE, method = "glm", aes(group=psTerminal, col=psTerminal)) + 
	theme_minimal() + scale_colour_manual(values = c("mediumaquamarine", "maroon4")) +
	theme(plot.title = element_text(family="Times", face="bold", size=20)) +
	ylab("Gray Matter Density") +
	labs(title = "Parietal Gray Matter Density (Male)", subtitle = mainstring)

stats_parietal_sa_M <- lmer(Parietal_SA ~ Age*psTerminal + (Age | bblid), df[df$sex == 2,])
coeffdf <- summary(stats_parietal_sa_M)$coefficients
estimates <- coeffdf[,"Estimate"]
stderrors <- coeffdf[,"Std. Error"]
tvalues <- coeffdf[,"t value"]
if (abs(tvalues[["Age"]]) > 2) { age_est <- paste0(round(estimates[["Age"]], digits=3), "*") } else { age_est <- round(estimates[["Age"]], digits=3) }
if (abs(tvalues[["psTerminalTD"]]) > 2) { diag_est <- paste0(round(estimates[["psTerminalTD"]], digits=3), "*") } else { diag_est <- round(estimates[["psTerminalTD"]], digits=3) }
if (abs(tvalues[["Age:psTerminalTD"]]) > 2) { inter_est <- paste0(round(estimates[["Age:psTerminalTD"]], digits=3), "*") } else { inter_est <- round(estimates[["Age:psTerminalTD"]], digits=3) }
mainstring <- paste0(malestring, "\n", "Coefficients: Age=", age_est, ", Diagnosis=", diag_est, ", Age:Diagnosis=", inter_est)
parietal_sa_M <- ggplot(df[df$sex == 2,], aes(x=Age, y=Parietal_SA)) + 
	geom_line(aes(group=bblid, col=psTerminal), alpha=.2) + 
	geom_smooth(se = TRUE, method = "glm", aes(group=psTerminal, col=psTerminal)) + 
	theme_minimal() + scale_colour_manual(values = c("mediumaquamarine", "maroon4")) +
	theme(plot.title = element_text(family="Times", face="bold", size=20)) +
	ylab("Surface Area (mm2)") +
	labs(title = "Parietal Surface Area (Male)", subtitle = mainstring)

# Temporal
stats_temporal_cort_M <- lmer(Temporal_CT ~ Age*psTerminal + (Age | bblid), df[df$sex == 2,])
coeffdf <- summary(stats_temporal_cort_M)$coefficients
estimates <- coeffdf[,"Estimate"]
stderrors <- coeffdf[,"Std. Error"]
tvalues <- coeffdf[,"t value"]
if (abs(tvalues[["Age"]]) > 2) { age_est <- paste0(round(estimates[["Age"]], digits=3), "*") } else { age_est <- round(estimates[["Age"]], digits=3) }
if (abs(tvalues[["psTerminalTD"]]) > 2) { diag_est <- paste0(round(estimates[["psTerminalTD"]], digits=3), "*") } else { diag_est <- round(estimates[["psTerminalTD"]], digits=3) }
if (abs(tvalues[["Age:psTerminalTD"]]) > 2) { inter_est <- paste0(round(estimates[["Age:psTerminalTD"]], digits=3), "*") } else { inter_est <- round(estimates[["Age:psTerminalTD"]], digits=3) }
mainstring <- paste0(malestring, "\n", "Coefficients: Age=", age_est, ", Diagnosis=", diag_est, ", Age:Diagnosis=", inter_est)
temporal_cort_M <- ggplot(df[df$sex == 2,], aes(x=Age, y=Temporal_CT)) + 
	geom_line(aes(group=bblid, col=psTerminal), alpha=.2) + 
	geom_smooth(se = TRUE, method = "glm", aes(group=psTerminal, col=psTerminal)) + 
	theme_minimal() + scale_colour_manual(values = c("mediumaquamarine", "maroon4")) +
	theme(plot.title = element_text(family="Times", face="bold", size=20)) +
	ylab("Cortical Thickness (mm)") +
	labs(title = "Temporal Cortical Thickness (Male)", subtitle = mainstring)

stats_temporal_vol_M <- lmer(Temporal_Vol ~ Age*psTerminal + (Age | bblid), df[df$sex == 2,])
coeffdf <- summary(stats_temporal_vol_M)$coefficients
estimates <- coeffdf[,"Estimate"]
stderrors <- coeffdf[,"Std. Error"]
tvalues <- coeffdf[,"t value"]
if (abs(tvalues[["Age"]]) > 2) { age_est <- paste0(round(estimates[["Age"]], digits=3), "*") } else { age_est <- round(estimates[["Age"]], digits=3) }
if (abs(tvalues[["psTerminalTD"]]) > 2) { diag_est <- paste0(round(estimates[["psTerminalTD"]], digits=3), "*") } else { diag_est <- round(estimates[["psTerminalTD"]], digits=3) }
if (abs(tvalues[["Age:psTerminalTD"]]) > 2) { inter_est <- paste0(round(estimates[["Age:psTerminalTD"]], digits=3), "*") } else { inter_est <- round(estimates[["Age:psTerminalTD"]], digits=3) }
mainstring <- paste0(malestring, "\n", "Coefficients: Age=", age_est, ", Diagnosis=", diag_est, ", Age:Diagnosis=", inter_est)
temporal_vol_M <- ggplot(df[df$sex == 2,], aes(x=Age, y=Temporal_Vol)) + 
	geom_line(aes(group=bblid, col=psTerminal), alpha=.2) + 
	geom_smooth(se = TRUE, method = "glm", aes(group=psTerminal, col=psTerminal)) + 
	theme_minimal() + scale_colour_manual(values = c("mediumaquamarine", "maroon4")) +
	theme(plot.title = element_text(family="Times", face="bold", size=20)) +
	ylab("Gray Matter Volume (mm3)") +
	labs(title = "Temporal Gray Matter Volume (Male)", subtitle = mainstring)

stats_temporal_gmd_M <- lmer(Temporal_GMD ~ Age*psTerminal + (Age | bblid), df[df$sex == 2,])
coeffdf <- summary(stats_temporal_gmd_M)$coefficients
estimates <- coeffdf[,"Estimate"]
stderrors <- coeffdf[,"Std. Error"]
tvalues <- coeffdf[,"t value"]
if (abs(tvalues[["Age"]]) > 2) { age_est <- paste0(round(estimates[["Age"]], digits=3), "*") } else { age_est <- round(estimates[["Age"]], digits=3) }
if (abs(tvalues[["psTerminalTD"]]) > 2) { diag_est <- paste0(round(estimates[["psTerminalTD"]], digits=3), "*") } else { diag_est <- round(estimates[["psTerminalTD"]], digits=3) }
if (abs(tvalues[["Age:psTerminalTD"]]) > 2) { inter_est <- paste0(round(estimates[["Age:psTerminalTD"]], digits=3), "*") } else { inter_est <- round(estimates[["Age:psTerminalTD"]], digits=3) }
mainstring <- paste0(malestring, "\n", "Coefficients: Age=", age_est, ", Diagnosis=", diag_est, ", Age:Diagnosis=", inter_est)
temporal_gmd_M <- ggplot(df[df$sex == 2,], aes(x=Age, y=Temporal_GMD)) + 
	geom_line(aes(group=bblid, col=psTerminal), alpha=.2) + 
	geom_smooth(se = TRUE, method = "glm", aes(group=psTerminal, col=psTerminal)) + 
	theme_minimal() + scale_colour_manual(values = c("mediumaquamarine", "maroon4")) +
	theme(plot.title = element_text(family="Times", face="bold", size=20)) +
	ylab("Gray Matter Density") +
	labs(title = "Temporal Gray Matter Density (Male)", subtitle = mainstring)

stats_temporal_sa_M <- lmer(Temporal_SA ~ Age*psTerminal + (Age | bblid), df[df$sex == 2,])
coeffdf <- summary(stats_temporal_sa_M)$coefficients
estimates <- coeffdf[,"Estimate"]
stderrors <- coeffdf[,"Std. Error"]
tvalues <- coeffdf[,"t value"]
if (abs(tvalues[["Age"]]) > 2) { age_est <- paste0(round(estimates[["Age"]], digits=3), "*") } else { age_est <- round(estimates[["Age"]], digits=3) }
if (abs(tvalues[["psTerminalTD"]]) > 2) { diag_est <- paste0(round(estimates[["psTerminalTD"]], digits=3), "*") } else { diag_est <- round(estimates[["psTerminalTD"]], digits=3) }
if (abs(tvalues[["Age:psTerminalTD"]]) > 2) { inter_est <- paste0(round(estimates[["Age:psTerminalTD"]], digits=3), "*") } else { inter_est <- round(estimates[["Age:psTerminalTD"]], digits=3) }
mainstring <- paste0(malestring, "\n", "Coefficients: Age=", age_est, ", Diagnosis=", diag_est, ", Age:Diagnosis=", inter_est)
temporal_sa_M <- ggplot(df[df$sex == 2,], aes(x=Age, y=Temporal_SA)) + 
	geom_line(aes(group=bblid, col=psTerminal), alpha=.2) + 
	geom_smooth(se = TRUE, method = "glm", aes(group=psTerminal, col=psTerminal)) + 
	theme_minimal() + scale_colour_manual(values = c("mediumaquamarine", "maroon4")) +
	theme(plot.title = element_text(family="Times", face="bold", size=20)) +
	ylab("Surface Area (mm2)") +
	labs(title = "Temporal Surface Area (Male)", subtitle = mainstring)

# Occipital
stats_occipital_cort_M <- lmer(Occipital_CT ~ Age*psTerminal + (Age | bblid), df[df$sex == 2,])
coeffdf <- summary(stats_occipital_cort_M)$coefficients
estimates <- coeffdf[,"Estimate"]
stderrors <- coeffdf[,"Std. Error"]
tvalues <- coeffdf[,"t value"]
if (abs(tvalues[["Age"]]) > 2) { age_est <- paste0(round(estimates[["Age"]], digits=3), "*") } else { age_est <- round(estimates[["Age"]], digits=3) }
if (abs(tvalues[["psTerminalTD"]]) > 2) { diag_est <- paste0(round(estimates[["psTerminalTD"]], digits=3), "*") } else { diag_est <- round(estimates[["psTerminalTD"]], digits=3) }
if (abs(tvalues[["Age:psTerminalTD"]]) > 2) { inter_est <- paste0(round(estimates[["Age:psTerminalTD"]], digits=3), "*") } else { inter_est <- round(estimates[["Age:psTerminalTD"]], digits=3) }
mainstring <- paste0(malestring, "\n", "Coefficients: Age=", age_est, ", Diagnosis=", diag_est, ", Age:Diagnosis=", inter_est)
occipital_cort_M <- ggplot(df[df$sex == 2,], aes(x=Age, y=Occipital_CT)) + 
	geom_line(aes(group=bblid, col=psTerminal), alpha=.2) + 
	geom_smooth(se = TRUE, method = "glm", aes(group=psTerminal, col=psTerminal)) + 
	theme_minimal() + scale_colour_manual(values = c("mediumaquamarine", "maroon4")) +
	theme(plot.title = element_text(family="Times", face="bold", size=20)) +
	ylab("Cortical Thickness (mm)") +
	labs(title = "Occipital Cortical Thickness (Male)", subtitle = mainstring)

stats_occipital_vol_M <- lmer(Occipital_Vol ~ Age*psTerminal + (Age | bblid), df[df$sex == 2,])
coeffdf <- summary(stats_occipital_vol_M)$coefficients
estimates <- coeffdf[,"Estimate"]
stderrors <- coeffdf[,"Std. Error"]
tvalues <- coeffdf[,"t value"]
if (abs(tvalues[["Age"]]) > 2) { age_est <- paste0(round(estimates[["Age"]], digits=3), "*") } else { age_est <- round(estimates[["Age"]], digits=3) }
if (abs(tvalues[["psTerminalTD"]]) > 2) { diag_est <- paste0(round(estimates[["psTerminalTD"]], digits=3), "*") } else { diag_est <- round(estimates[["psTerminalTD"]], digits=3) }
if (abs(tvalues[["Age:psTerminalTD"]]) > 2) { inter_est <- paste0(round(estimates[["Age:psTerminalTD"]], digits=3), "*") } else { inter_est <- round(estimates[["Age:psTerminalTD"]], digits=3) }
mainstring <- paste0(malestring, "\n", "Coefficients: Age=", age_est, ", Diagnosis=", diag_est, ", Age:Diagnosis=", inter_est)
occipital_vol_M <- ggplot(df[df$sex == 2,], aes(x=Age, y=Occipital_Vol)) + 
	geom_line(aes(group=bblid, col=psTerminal), alpha=.2) + 
	geom_smooth(se = TRUE, method = "glm", aes(group=psTerminal, col=psTerminal)) + 
	theme_minimal() + scale_colour_manual(values = c("mediumaquamarine", "maroon4")) +
	theme(plot.title = element_text(family="Times", face="bold", size=20)) +
	ylab("Gray Matter Volume (mm3)") +
	labs(title = "Occipital Gray Matter Volume (Male)", subtitle = mainstring)

stats_occipital_gmd_M <- lmer(Occipital_GMD ~ Age*psTerminal + (Age | bblid), df[df$sex == 2,])
coeffdf <- summary(stats_occipital_gmd_M)$coefficients
estimates <- coeffdf[,"Estimate"]
stderrors <- coeffdf[,"Std. Error"]
tvalues <- coeffdf[,"t value"]
if (abs(tvalues[["Age"]]) > 2) { age_est <- paste0(round(estimates[["Age"]], digits=3), "*") } else { age_est <- round(estimates[["Age"]], digits=3) }
if (abs(tvalues[["psTerminalTD"]]) > 2) { diag_est <- paste0(round(estimates[["psTerminalTD"]], digits=3), "*") } else { diag_est <- round(estimates[["psTerminalTD"]], digits=3) }
if (abs(tvalues[["Age:psTerminalTD"]]) > 2) { inter_est <- paste0(round(estimates[["Age:psTerminalTD"]], digits=3), "*") } else { inter_est <- round(estimates[["Age:psTerminalTD"]], digits=3) }
mainstring <- paste0(malestring, "\n", "Coefficients: Age=", age_est, ", Diagnosis=", diag_est, ", Age:Diagnosis=", inter_est)
occipital_gmd_M <- ggplot(df[df$sex == 2,], aes(x=Age, y=Occipital_GMD)) + 
	geom_line(aes(group=bblid, col=psTerminal), alpha=.2) + 
	geom_smooth(se = TRUE, method = "glm", aes(group=psTerminal, col=psTerminal)) + 
	theme_minimal() + scale_colour_manual(values = c("mediumaquamarine", "maroon4")) +
	theme(plot.title = element_text(family="Times", face="bold", size=20)) +
	ylab("Gray Matter Density") +
	labs(title = "Occipital Gray Matter Density (Male)", subtitle = mainstring)

stats_occipital_sa_M <- lmer(Occipital_SA ~ Age*psTerminal + (Age | bblid), df[df$sex == 2,])
coeffdf <- summary(stats_occipital_sa_M)$coefficients
estimates <- coeffdf[,"Estimate"]
stderrors <- coeffdf[,"Std. Error"]
tvalues <- coeffdf[,"t value"]
if (abs(tvalues[["Age"]]) > 2) { age_est <- paste0(round(estimates[["Age"]], digits=3), "*") } else { age_est <- round(estimates[["Age"]], digits=3) }
if (abs(tvalues[["psTerminalTD"]]) > 2) { diag_est <- paste0(round(estimates[["psTerminalTD"]], digits=3), "*") } else { diag_est <- round(estimates[["psTerminalTD"]], digits=3) }
if (abs(tvalues[["Age:psTerminalTD"]]) > 2) { inter_est <- paste0(round(estimates[["Age:psTerminalTD"]], digits=3), "*") } else { inter_est <- round(estimates[["Age:psTerminalTD"]], digits=3) }
mainstring <- paste0(malestring, "\n", "Coefficients: Age=", age_est, ", Diagnosis=", diag_est, ", Age:Diagnosis=", inter_est)
occipital_sa_M <- ggplot(df[df$sex == 2,], aes(x=Age, y=Occipital_SA)) + 
	geom_line(aes(group=bblid, col=psTerminal), alpha=.2) + 
	geom_smooth(se = TRUE, method = "glm", aes(group=psTerminal, col=psTerminal)) + 
	theme_minimal() + scale_colour_manual(values = c("mediumaquamarine", "maroon4")) +
	theme(plot.title = element_text(family="Times", face="bold", size=20)) +
	ylab("Surface Area (mm2)") +
	labs(title = "Occipital Surface Area (Male)", subtitle = mainstring)

# BasGang
stats_basGang_vol_M <- lmer(BasGang_Vol ~ Age*psTerminal + (Age | bblid), df[df$sex == 2,])
coeffdf <- summary(stats_basGang_vol_M)$coefficients
estimates <- coeffdf[,"Estimate"]
stderrors <- coeffdf[,"Std. Error"]
tvalues <- coeffdf[,"t value"]
if (abs(tvalues[["Age"]]) > 2) { age_est <- paste0(round(estimates[["Age"]], digits=3), "*") } else { age_est <- round(estimates[["Age"]], digits=3) }
if (abs(tvalues[["psTerminalTD"]]) > 2) { diag_est <- paste0(round(estimates[["psTerminalTD"]], digits=3), "*") } else { diag_est <- round(estimates[["psTerminalTD"]], digits=3) }
if (abs(tvalues[["Age:psTerminalTD"]]) > 2) { inter_est <- paste0(round(estimates[["Age:psTerminalTD"]], digits=3), "*") } else { inter_est <- round(estimates[["Age:psTerminalTD"]], digits=3) }
mainstring <- paste0(malestring, "\n", "Coefficients: Age=", age_est, ", Diagnosis=", diag_est, ", Age:Diagnosis=", inter_est)
basGang_vol_M <- ggplot(df[df$sex == 2,], aes(x=Age, y=BasGang_Vol)) + 
	geom_line(aes(group=bblid, col=psTerminal), alpha=.2) + 
	geom_smooth(se = TRUE, method = "glm", aes(group=psTerminal, col=psTerminal)) + 
	theme_minimal() + scale_colour_manual(values = c("mediumaquamarine", "maroon4")) +
	theme(plot.title = element_text(family="Times", face="bold", size=20)) +
	ylab("Gray Matter Volume (mm3)") +
	labs(title = "Basal Gangalia Gray Matter Volume (Male)", subtitle = mainstring)

stats_basGang_gmd_M <- lmer(BasGang_GMD ~ Age*psTerminal + (Age | bblid), df[df$sex == 2,])
coeffdf <- summary(stats_basGang_gmd_M)$coefficients
estimates <- coeffdf[,"Estimate"]
stderrors <- coeffdf[,"Std. Error"]
tvalues <- coeffdf[,"t value"]
if (abs(tvalues[["Age"]]) > 2) { age_est <- paste0(round(estimates[["Age"]], digits=3), "*") } else { age_est <- round(estimates[["Age"]], digits=3) }
if (abs(tvalues[["psTerminalTD"]]) > 2) { diag_est <- paste0(round(estimates[["psTerminalTD"]], digits=3), "*") } else { diag_est <- round(estimates[["psTerminalTD"]], digits=3) }
if (abs(tvalues[["Age:psTerminalTD"]]) > 2) { inter_est <- paste0(round(estimates[["Age:psTerminalTD"]], digits=3), "*") } else { inter_est <- round(estimates[["Age:psTerminalTD"]], digits=3) }
mainstring <- paste0(malestring, "\n", "Coefficients: Age=", age_est, ", Diagnosis=", diag_est, ", Age:Diagnosis=", inter_est)
basGang_gmd_M <- ggplot(df[df$sex == 2,], aes(x=Age, y=BasGang_GMD)) + 
	geom_line(aes(group=bblid, col=psTerminal), alpha=.2) + 
	geom_smooth(se = TRUE, method = "glm", aes(group=psTerminal, col=psTerminal)) + 
	theme_minimal() + scale_colour_manual(values = c("mediumaquamarine", "maroon4")) +
	theme(plot.title = element_text(family="Times", face="bold", size=20)) +
	ylab("Gray Matter Density") +
	labs(title = "Basal Gangalia Gray Matter Density (Male)", subtitle = mainstring)

# Limbic
stats_limbic_vol_M <- lmer(Limbic_Vol ~ Age*psTerminal + (Age | bblid), df[df$sex == 2,])
coeffdf <- summary(stats_limbic_vol_M)$coefficients
estimates <- coeffdf[,"Estimate"]
stderrors <- coeffdf[,"Std. Error"]
tvalues <- coeffdf[,"t value"]
if (abs(tvalues[["Age"]]) > 2) { age_est <- paste0(round(estimates[["Age"]], digits=3), "*") } else { age_est <- round(estimates[["Age"]], digits=3) }
if (abs(tvalues[["psTerminalTD"]]) > 2) { diag_est <- paste0(round(estimates[["psTerminalTD"]], digits=3), "*") } else { diag_est <- round(estimates[["psTerminalTD"]], digits=3) }
if (abs(tvalues[["Age:psTerminalTD"]]) > 2) { inter_est <- paste0(round(estimates[["Age:psTerminalTD"]], digits=3), "*") } else { inter_est <- round(estimates[["Age:psTerminalTD"]], digits=3) }
mainstring <- paste0(malestring, "\n", "Coefficients: Age=", age_est, ", Diagnosis=", diag_est, ", Age:Diagnosis=", inter_est)
limbic_vol_M <- ggplot(df[df$sex == 2,], aes(x=Age, y=Limbic_Vol)) + 
	geom_line(aes(group=bblid, col=psTerminal), alpha=.2) + 
	geom_smooth(se = TRUE, method = "glm", aes(group=psTerminal, col=psTerminal)) + 
	theme_minimal() + scale_colour_manual(values = c("mediumaquamarine", "maroon4")) +
	theme(plot.title = element_text(family="Times", face="bold", size=20)) +
	ylab("Gray Matter Volume (mm3)") +
	labs(title = "Limbic Gray Matter Volume (Male)", subtitle = mainstring)

stats_limbic_gmd_M <- lmer(Limbic_GMD ~ Age*psTerminal + (Age | bblid), df[df$sex == 2,])
coeffdf <- summary(stats_limbic_gmd_M)$coefficients
estimates <- coeffdf[,"Estimate"]
stderrors <- coeffdf[,"Std. Error"]
tvalues <- coeffdf[,"t value"]
if (abs(tvalues[["Age"]]) > 2) { age_est <- paste0(round(estimates[["Age"]], digits=3), "*") } else { age_est <- round(estimates[["Age"]], digits=3) }
if (abs(tvalues[["psTerminalTD"]]) > 2) { diag_est <- paste0(round(estimates[["psTerminalTD"]], digits=3), "*") } else { diag_est <- round(estimates[["psTerminalTD"]], digits=3) }
if (abs(tvalues[["Age:psTerminalTD"]]) > 2) { inter_est <- paste0(round(estimates[["Age:psTerminalTD"]], digits=3), "*") } else { inter_est <- round(estimates[["Age:psTerminalTD"]], digits=3) }
mainstring <- paste0(malestring, "\n", "Coefficients: Age=", age_est, ", Diagnosis=", diag_est, ", Age:Diagnosis=", inter_est)
limbic_gmd_M <- ggplot(df[df$sex == 2,], aes(x=Age, y=Limbic_GMD)) + 
	geom_line(aes(group=bblid, col=psTerminal), alpha=.2) + 
	geom_smooth(se = TRUE, method = "glm", aes(group=psTerminal, col=psTerminal)) + 
	theme_minimal() + scale_colour_manual(values = c("mediumaquamarine", "maroon4")) +
	theme(plot.title = element_text(family="Times", face="bold", size=20)) +
	ylab("Gray Matter Density") +
	labs(title = "Limbic Gray Matter Density (Male)", subtitle = mainstring)



# export plots
pdf(file="/home/butellyn/longitudinal_psychosis/plots/F_lobular_longitudinal_lmer_psTerminal.pdf", width=14, height=6)
grid.arrange(frontorb_cort_F, frontorb_vol_F, ncol=2)
grid.arrange(frontorb_gmd_F, front_sa_F, ncol=2)
grid.arrange(frontdors_cort_F, frontdors_vol_F, ncol=2)
grid.arrange(frontdors_gmd_F, front_sa_F, ncol=2)
grid.arrange(parietal_cort_F, parietal_vol_F, ncol=2)
grid.arrange(parietal_gmd_F, parietal_sa_F, ncol=2)
grid.arrange(temporal_cort_F, temporal_vol_F, ncol=2)
grid.arrange(temporal_gmd_F, temporal_sa_F, ncol=2)
grid.arrange(occipital_cort_F, occipital_vol_F, ncol=2)
grid.arrange(occipital_gmd_F, occipital_sa_F, ncol=2)
grid.arrange(limbic_vol_F, limbic_gmd_F, ncol=2)
grid.arrange(basGang_vol_F, basGang_gmd_F, ncol=2)
dev.off()

pdf(file="/home/butellyn/longitudinal_psychosis/plots/M_lobular_longitudinal_lmer_psTerminal.pdf", width=14, height=6)
grid.arrange(frontorb_cort_M, frontorb_vol_M, ncol=2)
grid.arrange(frontorb_gmd_M, front_sa_M, ncol=2)
grid.arrange(frontdors_cort_M, frontdors_vol_M, ncol=2)
grid.arrange(frontdors_gmd_M, front_sa_M, ncol=2)
grid.arrange(parietal_cort_M, parietal_vol_M, ncol=2)
grid.arrange(parietal_gmd_M, parietal_sa_M, ncol=2)
grid.arrange(temporal_cort_M, temporal_vol_M, ncol=2)
grid.arrange(temporal_gmd_M, temporal_sa_M, ncol=2)
grid.arrange(occipital_cort_M, occipital_vol_M, ncol=2)
grid.arrange(occipital_gmd_M, occipital_sa_M, ncol=2)
grid.arrange(limbic_vol_M, limbic_gmd_M, ncol=2)
grid.arrange(basGang_vol_M, basGang_gmd_M, ncol=2)
dev.off()











