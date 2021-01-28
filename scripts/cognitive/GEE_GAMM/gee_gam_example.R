# Example using GAM + GEE to test significance on example vertices and in simulated data

library(mgcv)
library(geepack)
# install.packages("doBy")
library(doBy)
library(MASS)

setwd("/home/smweinst/test_vertices_gee_gam/")

## Example vertices ----
# load in example vertices
Long_Motor<-readRDS('LongFormat_MotorVert.rds')
Long_Visual<-readRDS('LongFormat_VisualVert.rds')
Long_PFC<-readRDS('LongFormat_PFCVert.rds')
# load in mean vertex values at each scale for each subj. as pseudo example vertex
Long_Avg<-readRDS('LongFormat_AvgVert.rds')

# is "value" column of Long_PFC the same as BW_FC in the other datasets? going to assume yes (need consistent variable names for the loop below)
names(Long_PFC)[which(names(Long_PFC)=="value")] = "BW_FC"

# list different vertex datasets
list_vertices = list(motor = Long_Motor, visual = Long_Visual, pfc = Long_PFC, avg = Long_Avg)
n_vertex = length(list_vertices)

# L_contrast:
# 0's correspond to coefficients where variable is included in both the full and reduced model
# 1's correspond to coefficients that we want to test (jointly) if = 0
L_contrast = c(rep(0,3),rep(1,8))

pvalX2 = vector(mode = "numeric", length = n_vertex)

# p-value for every vertex
names(pvalX2) =  names(list_vertices)

for (v in 1:n_vertex){ # loop over the vertices
  vertex_dat = list_vertices[[v]]
  
  # gam with smooth terms for scale, age, and interaction between scale and age
  # doesn't account for within-subject correlation, but running this model so that we can extract the model.matrix
  fit_gam = mgcv::gam(BW_FC~Motion+Sex+s(Scale,k=3,fx=T)+s(Age,k=3,fx=T)+ti(Scale,Age,k=3,fx=T),
                      data=vertex_dat)
  
  # extract model matrix from GAM and then input the smooth terms to a gee (assuming exchangeable correlation structure, might be wrong)
  gam.model.matrix = cbind(model.matrix(fit_gam), bblid = vertex_dat$bblid,
                           BW_FC = vertex_dat$BW_FC, Scale = vertex_dat$Scale)
  gam.model.df = data.frame(gam.model.matrix)
  gam.model.df = gam.model.df[order(gam.model.df$bblid),] # sort by ID for geeglm
  
  # doing the rest with a geeglm will give robust variance estimators, accounting for within-subject
  # correlation which was not the case in the gam part of the output from GAMM when we were doing bootstrapping
  fit_gee = geepack::geeglm(BW_FC~Motion + Sex2 + s.Scale..1 + 
                              s.Scale..2 + s.Age..1 + s.Age..2 + 
                              ti.Scale.Age..1 + ti.Scale.Age..2 + ti.Scale.Age..3 +
                              ti.Scale.Age..4, id = gam.model.df$bblid,
                            data = gam.model.df, corstr = "exchangeable")
  
  # joint test of whether the coefficients corresponding to a `1` in L_contrast differ from 0
  joint_test = doBy::esticon(obj = fit_gee,L = L_contrast,joint.test = T)
  
  # p-value based on chi-square with 1 df
  pvalX2[v] = pchisq(joint_test$X2.stat, df = 1, lower.tail = F)
  
}

# p-value for each example vertex:
# print(pvalX2)
#     motor        visual           pfc           avg 
# 3.615059e-53  1.487835e-17  4.321208e-05 9.509597e-184 

# example testing a single coefficient instead of doing a joint test:
vertex_dat = list_vertices[[3]]

fit_gam = mgcv::gam(BW_FC~Motion+Sex+s(Scale,k=3,fx=T)+s(Age,k=3,fx=T)+ti(Scale,Age,k=3,fx=T),
                    data=vertex_dat)
gam.model.matrix = cbind(model.matrix(fit_gam), bblid = vertex_dat$bblid,
                         BW_FC = vertex_dat$BW_FC)
gam.model.df = data.frame(gam.model.matrix)
gam.model.df = gam.model.df[order(gam.model.df$bblid),] # sort by ID for geeglm

# doing the rest with a geeglm will give robust variance estimators, accounting for within-subject
# correlation which was not the case in the gam part of the output from GAMM when we were doing bootstrapping
fit_gee = geepack::geeglm(BW_FC~Motion + Sex2 + s.Scale..1 + 
                            s.Scale..2 + s.Age..1 + s.Age..2 + 
                            ti.Scale.Age..1 + ti.Scale.Age..2 + ti.Scale.Age..3 +
                            ti.Scale.Age..4, id = as.factor(gam.model.df$bblid),
                          data = gam.model.df, corstr = "exchangeable")



# joint test of whether the coefficients corresponding to a `1` in L_contrast differ from 0
joint_test_TRUE = doBy::esticon(obj = fit_gee,L = c(0,1,rep(0,9)),joint.test = T)

print(joint_test_TRUE)

joint_test_FALSE = doBy::esticon(obj = fit_gee,L = c(0,1,rep(0,9)),joint.test = F)
print(joint_test_FALSE)

## Simulations ----

# X1 ~ N(0,1) for subjects 1,...,n (one per subject)

# (v1, v2, v3) ~ MVN((0,0,0), Sigma)     Sigma is a 3x3 covariance matrix
## transform data from wide to long format so that the column "v" in the long data includes v1, v2, and v3 for each subject
## (this is supposed to be analogous to the "scale" variable in the vertex examples)

# correlation structure scenarios:
## (1) correlation structure is exchangeable (correct specification of the GEE)
## and (2) correlation structure is not exchangeable (incorrect specification of the GEE)

# outcome variable:
# Y = a*(X1) + b*(v^3 + (v^3)*X1) + rnorm(n) -- add N(0,1) noise to each observation
# a and b are coefficients specified in the simulations; we want to test whether b = 0 (joint test of whether the coefficient(s) on the smooth terms from the GEE are 0)


# covariance matrices:
## exchangeable correlation/covariance structure:
cov_empty =  matrix(nrow=3,ncol=3)
diag(cov_empty) = 1
cov1 = cov_empty; cov1[which(is.na(cov1))] = 0.9
cov2 = cov_empty; cov2[which(is.na(cov2))] = 0.6

# one example with some other correlation structure (gee correlation misspecified):
cov3 = cov_empty; cov3[upper.tri(cov3)] = c(0.2,0.6,0.1); cov3[lower.tri(cov3)] = cov3[upper.tri(cov3)]

par.vals = expand.grid(a = 1, # coefficient for linear terms in true model for Y
                       b = seq(0,0.1,by = 0.02), # coefficient for smooth terms in true model for Y
                       # null hypothesis of interest is true when b = 0
                       covar = paste0("cov", 1:3)
)

nsim = 1000

set.seed(917)

out = list()
pb <- txtProgressBar(min = 0, max = nrow(par.vals), style = 3)

system.time(
for (j in 1:nrow(par.vals)){
  out[[j]] = vector(mode = "numeric", length = nsim)
  
  for (sim in 1:nsim){
    dat.sim = MASS::mvrnorm(n = n, mu = rep(0,3), Sigma = get(as.character(par.vals$covar[j])))
    colnames(dat.sim) = paste0("v",1:3)
    dat.sim = data.frame(cbind(id = 1:n, dat.sim))
    dat.sim$X1 = rnorm(n) # covariate X1
   # dat.sim$X2 = rnorm(n)
    
    # convert data from wide to long format (multiple rows per subject)
    dat.sim_long = reshape(dat.sim, idvar = "id", varying = list(c("v1","v2","v3")),v.names = c("v"),direction = "long")
    dat.sim_long = dat.sim_long[order(dat.sim_long$id),]
    
    dat.sim_long$Y = par.vals$a[j]*(dat.sim_long$X1 ) + par.vals$b[j]*(dat.sim_long$v^3 + (dat.sim_long$v^3)*dat.sim_long$X1) + rnorm(n) # rnorm(n) for additional noise
    
    # fit gam to get model.matrix
    fit_gam = mgcv::gam(Y~X1+s(v,k=3,fx=T)+ ti(X1,v,k=3,fx=T),
                        data=dat.sim_long)
    gam.model.matrix = cbind(model.matrix(fit_gam), id = dat.sim_long$id,
                             Y = dat.sim_long$Y)
    
    gam.model.df = data.frame(gam.model.matrix)
    
    # gee using terms from model.matrix from the gam above; assuming exchangeable correlation structure (which is true in the simulated data...)
    fit_gee = geepack::geeglm(Y~X1 + s.v..1 + s.v..2 + ti.X1.v..1 + ti.X1.v..2 + ti.X1.v..3 + ti.X1.v..4, id = gam.model.df$id,
                              data = gam.model.df, corstr = "exchangeable")

    # joint test of the smooth terms
    joint_test = doBy::esticon(obj = fit_gee,L = c(0,0,1,1,1,1,1,1),joint.test = T)
    out[[j]][sim] = pchisq(joint_test$X2.stat, df = 1, lower.tail = F)
    
  }
  
  par.vals$p_reject[j] = length(which(out[[j]]<0.05))/nsim # for each set of parameters, get the proportion of simulations where which H_0 is rejected
  setTxtProgressBar(pb, j)
  
}
)
close(pb)


# plot power and type I error
par(mfrow = c(1,1))
plot(x = par.vals$b, y=par.vals$p_reject, type = 'n',
     xlab = "b", ylab = "pr(reject joint H0)",ylim = c(0,1), xlim = c(min(par.vals$b), max(par.vals$b)))
abline(h = c(0.05,0.1), col = c("blue","red"), lwd = c(2,2))

# specs for plot by a and covariance
legend_specs = unique(par.vals[,c("a","covar")])
legend_specs$lty = as.numeric(substr(legend_specs$covar,4,4))

for (c in unique(par.vals$covar)){
    points(par.vals$b[which(par.vals$covar==c)],
           par.vals$p_reject[which(par.vals$covar==c)],
           type = 'b',
           lty = legend_specs$lty[which(legend_specs$covar==c)],
           lwd = 1.5, pch = 16, cex =0.8)
}


legend("topleft",lty = c(legend_specs$lty,1,1),
       legend = c("exchangeable (rho = 0.2)", "exchangeable (rho = 0.6)", "misspecified correlation structure", "p = 0.05", "p = 0.1"),
       col = c(rep("black",3),"blue","red"),lwd = 1.5, bty = 'n',
       cex = 0.8)

# cov1 - exchangeable
#       [,1] [,2] [,3]
# [1,]  1.0  0.2  0.2
# [2,]  0.2  1.0  0.2
# [3,]  0.2  0.2  1.0

# cov2 - exchangeable
#       [,1] [,2] [,3]
# [1,]  1.0  0.6  0.6
# [2,]  0.6  1.0  0.6
# [3,]  0.6  0.6  1.0

# cov3 - misspecified
#     [,1] [,2] [,3]
# [1,]  1.0  0.2  0.6
# [2,]  0.2  1.0  0.1
# [3,]  0.6  0.1  1.0

