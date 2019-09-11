library('R.matlab')
library('matrixcalc')
library("Matrix")
library("lme4")
setwd('C:\\Users\\yancx\\Desktop\\Thesis\\Peter fMRI data\\')
# basis <- readMat("basis.mat")
# load("Scans.arr")
#scan_25 <- Scans.arr[,,,1:25]
# save(scan_25, file = "scan_25")
# 
load("scan_25")
dim(scan_25)

plot(scan_25[,1,1,1], type = "l")
plot(scan_25[,2,1,1], type = "l")

cor(scan_25[,1,1,1],scan_25[,4,1,1])



library("arrayhelpers")
# a <- arrayhelpers:::a
# dim(a)

# array2df(a, matrix=T)

# test <- scan_25[1:5,1:4,1:3,1:3]

# test.mat <- array2df(test, matrix = T)
# colnames(test.mat) <- c("y","ts","reg","scan","ID")
# test.mat

# scan_25.mat <- array2df(scan_25, matrix = T)
# colnames(scan_25.mat) <- c("y","ts","reg","sca","ID")
# scan_25.df <- data.frame(scan_25.mat)

# save(scan_25.df, file = "scan_25.df")
load("scan_10.df")
head(scan_10.df)


fm1 <- lmer(y ~ reg + sca + ID + (1 + reg|ID), data = scan_10.df)


summary(fm1)

t1 <- expand.grid(x=1:116, y= 1:116)


