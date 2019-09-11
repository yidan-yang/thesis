library('R.matlab')
library('matrixcalc')
setwd('C:\\Users\\yancx\\Desktop\\Thesis\\Peter fMRI data\\')
basis <- readMat("basis.mat")
load("Scans.arr")


dim(Scans.arr)

# One scan
Y <- Scans.arr[,,1,]

# One scan one patient
Y1 <- Y[,,1]

dimreformat2 <- function(Y){
  h <- 100
  k <- dim(Y)[2]
  t <- dim(Y)[1]-h+1
  # aa1 is for one patient MA time series
  # aa1 will eventually convert into matrix then store in an array aa2
  # finally it will convert into matrix w/ dim(h,t*k)
  aa1 <- array(data = 0,dim = c(h,k,t))
  # aa2 is for all 820 patients
  # also the final output array w/ dim 100*(1101*116)*820
  aa2 <- array(data = 0,dim = c(h,t*k,dim(Y)[3]))
  for (ind in 1:dim(Y)[3]) {
    for (i in 1:t) {
      aa1[,,i] <- Y[i:(i+99),,ind]
      aa2[,,ind] <- matrix(aa1,ncol = t*k)
    }
  }
  aa2
}
Y_30s <-dimreformat2(Y)
save("Y_30s", file = Y_30s)


#######################################
library(fda)
time_length=30
Breaks=0:time_length
ttime <- 1:time_length

B<-bsplineS(ttime,breaks=Breaks,norder=4)
r <- 7

D<-matrix(0,nrow=time_length,ncol=r)

for( r in 1:7){
  D[,r] <- sqrt(2/time_length)*cos(r*pi/time_length*(1:time_length))
}

matplot(B, type = "l", main = "B")
matplot(D, type = "l", main = "D")