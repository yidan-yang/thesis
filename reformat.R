library('R.matlab')
library('matrixcalc')
setwd('C:\\Users\\yancx\\Desktop\\Thesis\\Peter fMRI data\\')
basis <- readMat("basis.mat")
load("Scans.arr")
scan <- Scans.arr


Y <- scan[,,1,]

dim(Y)
dim(Y)[2]


# For one patient Y[,,1]
reformat <- function(Y){
  h <- 30
  k <- dim(Y)[2]
  t <- dim(Y)[1]-h+1
  aa1 <- array(data = 0,dim = c(h,k,t))
  Yn <- Y[,,1]
  for (i in 1:t) {
    aa1[,,i] <- Yn[i:(i+29),]
  }
  matrix(aa1, ncol = k*t)
}
reformat(Y)

dim(reformat(Y))
class(reformat(Y))


dim(Y)[3]

# For all 860 patients
reformat2 <- function(Y){
  h <- 30
  k <- dim(Y)[2]
  t <- dim(Y)[1]-h+1
  # aa1 is for one patient
  # finally it will convert into matrix w/ dim(h,t*k)
  aa1 <- array(data = 0,dim = c(h,k,t))
  # aa2 is for all 820 patients
  aa2 <- array(data = 0,dim = c(h,t*k,dim(Y)[3]))
  for (ind in 1:dim(Y)[3]) {
    Yn <- Y[,,ind]
    for (i in 1:t) {
      aa1[,,i] <- Yn[i:(i+29),]
    }
    matrix(aa1, ncol = k*t)
  }
}


#####################################
aa1 <- array(data = 0,dim = c(2,2,5))
mm2 <- matrix(data = rep(1:12), ncol = 2, byrow = F)

aa1[,,1] <- mm2[1:2,]
aa1[,,2] <- mm2[2:(2+1),]
aa1[,,3] <- mm2[3:4,]
aa1[,,4] <- mm2[4:5,]
aa1[,,5] <- mm2[5:6,]
aa1

aa1 <- array(data = 0,dim = c(2,2,5))
mm2 <- matrix(data = rep(1:12), ncol = 2, byrow = F)

for (i in 1:5) {
  aa1[,,i] <- mm2[i:(i+1),]
}
aa1

mtr <- function(mm2){
  aa1 <- array(data = 0,dim = c(2,2,5))
  mm2 <- matrix(data = rep(1:12), ncol = 2, byrow = F)
  for (i in 1:5) {
      aa1[,,i] <- mm2[i:(i+1),]
    }
  aa1
}
mtr(mm2)


######################## h*(t*k) matrix



aa1 <- array(data = 0,dim = c(2,2,5))
mm2 <- matrix(data = rep(1:12), ncol = 2, byrow = F)

for (i in 1:5) {
  aa1[,,i] <- mm2[i:(i+1),]
}
aa1

matrix(aa1, ncol = 10)


######################## multiple matrix store in an array test 1

# aa1 is for MA time series split
# aa1 will eventually convert into matrix then store in an array ma1
aa1 <- array(data = 0,dim = c(2,2,5))
aa1

# ta1 is similar to Y with 3D
ta1 <- array(data = rep(1:36), dim = c(6,2,3))
ta1

# mat1 is the final output array w/ dim 2*10 *3
ma1 <- array(data = 0, dim = c(2,10,3))
ma1

# in this example, h=2, t=6-2+1=5,k=2

# for 1 in ta1
for (i in 1:5) {
  aa1[,,i] <- ta1[i:(i+1),,1]
}
aa1
ma1[,,1]<- matrix(aa1, ncol = 10)
ma1


aa1 <- array(data = 0,dim = c(2,2,5))
ta1 <- array(data = rep(1:36), dim = c(6,2,3))
ma1 <- array(data = 0, dim = c(2,10,3))
# for all 3 in ta1
for (ind in 1:3) {
  for (i in 1:5) {
    aa1[,,i] <- ta1[i:(i+1),,ind]
    ma1[,,ind] <- matrix(aa1,ncol = 10)
  }
}
ma1

