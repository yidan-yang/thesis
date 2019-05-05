library('R.matlab')
library('matrixcalc')
setwd('C:\\Users\\yancx\\Desktop\\Thesis\\Peter fMRI data\\')
basis <- readMat("basis.mat")
load("Scans.arr")

scan <- Scans.arr

Y <- scan[,,1,]
dim(Y)
# Y is a list with 820 lists and each list is a 1200*116 matrix

head(Y)
dim(scan)
head(scan)

str(Y)

mat1 <- matrix(data = c(2,-2,1,-1,3,-1,2,-4,1),ncol = 3,byrow = T)
mat1
frobenius.norm(mat1)

Y

m1 <- matrix(data = c(1,2,3,4),ncol = 2, byrow = T)
m1
m2 <- matrix(data = c(2,0,1,2),ncol = 2, byrow = T)
m2
m1%*%m2

# d is random array with dimension (7,116,820)
# u is random matrix with dimension (33,p) where p=2 in this case
# v is random array with dimension (p,116,820) where p=2

p <- 2

d <- array(data = runif(7*116*820), dim = c(7,116,820))
u <- matrix(data = runif(33*p), nrow = 33, ncol = p)
v <- array(data = runif(p*116*820), dim = c(p,116,820))


test1 <- Y[,,1]-basis$D%*%d[,,1]-basis$B%*%u%*%v[,,1]
test1
class(test1)
dim(test1)
frobenius.norm(test1)
frobenius.norm(Y[,,1]-basis$D%*%d[,,1]-basis$B%*%u%*%v[,,1])

d1 <- array(data = 1:12, dim = c(2,3,2))
d1
mat11 <- matrix(data = c(1,1,1,0,0,0), ncol = 3, byrow = T)
mat11

out.mat1 <- NULL
func1 <- function(ind){
  for (ind in 1:2) {
    out = d1[,,ind] + mat11
    out.mat1 <- rbind(out.mat1,out)
  }
  out.mat1
}
func1(2)

frob.mat <- NULL
frob.func <- function(n){
  for (i in 1:n) {
    out =frobenius.norm(Y[,,i]-basis$D%*%d[,,i]-basis$B%*%u%*%v[,,i])   
    frob.mat <- rbind(frob.mat,out)
  }
  frob.mat
}
t11 <- frob.func(820)
class(t11)
dim(t11)
