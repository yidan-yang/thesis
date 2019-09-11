library('R.matlab')
library('matrixcalc')
library("Matrix")
library("lme4")
setwd('C:\\Users\\yancx\\Desktop\\Thesis\\Peter fMRI data\\')

load("scan_25")

# Since we have 116 regions, so we have C(116,2)=6670 combinations.
# And we create correlation matrix for 25 subjects, each subject has 4 scans.
# dim(cor.mat)=6670*(25*4)

reg.comb <- expand.grid(x=1:116,y=1:116)
reg.comb <- unique(t(apply(reg.comb,1,sort)))
reg.comb <- reg.comb[reg.comb[,1]!=reg.comb[,2],,drop = F]

cor.mat <- matrix(data = NA, nrow = 6670, ncol = 102)
cor.mat[,1:2] <- reg.comb

# Calculating correlation
dim(scan_25)
