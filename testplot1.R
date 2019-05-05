library('R.matlab')
library('matrixcalc')
setwd('C:\\Users\\yancx\\Desktop\\Thesis\\Peter fMRI data\\')
basis <- readMat("basis.mat")
load("Scans.arr")
scan <- Scans.arr

dim(scan)


Y <- scan[,,1,]


Yone <- Y[,,1]
dim(Yone)

win.graph()
for (i in 1:dim(Yone)[2]) {
  par(mfrow=c(dim(Yone)[2],1))
  ts.plot(Yone[,i])
}
