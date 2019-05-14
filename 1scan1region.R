library('R.matlab')
library('matrixcalc')
setwd('C:\\Users\\yancx\\Desktop\\Thesis\\Peter fMRI data\\')
basis <- readMat("basis.mat")
load("Scans.arr")
scan <- Scans.arr

dim(scan)


Y <- scan[,,1,]

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


R <- basis$R

dim(D)

p <- 2
v <- array(data = runif(p*1171), dim = c(p,1171,1))
d <- array(data = runif(7*1171), dim = c(7,1171,1))
u <- matrix(data = runif(33*p), nrow = 33, ncol = p)

# For patient 1, 1 scan and 1 region
Y1 <- Y[,,1]
Y_11 <- Y1[,1]
Y_11 <- as.matrix(Y_11)

# reposition the data to become 30*1171

reformat <- function(Y){
  h <- 30
  t <- dim(Y)[1]-h+1
  mat1 <- matrix(data = 0, nrow = h, ncol = t)
  Yn <- Y
  for (i in 1:t) {
    mat1[,i] <- Yn[i:(i+29),]
  }
  mat1
}
one_SandR <- reformat(Y_11)

one_arr <- array(data = one_SandR, dim = c(30,1171,1))


###################################################
#       Obj.val
###################################################
obj_val <- function(Y,D,d, B,u,v, lambda, R){
  n = dim(Y)[3]
  val=0;
  for (i in 1:n) {
    val = val + frobenius.norm(Y[,,i]-D%*%d[,,i]-B%*%u%*%v[,,i])^2   
  }
  frob.val = val/n +lambda*matrix.trace(t(u)%*%(R%*%u))
  frob.val
}

###################################################

# step1: fix U and d, solve v(v is unknown here)
# step2: fix v(at solution of previous step), u and d are unknown here

library("Matrix") # sparse matrix
library("pcg")
library("matlib")

lr.func.stp12 <- function(Y, D, d, B, u, v, lambda, R){
  # [T,L]=size(X)
  mat.TL <- matrix(c(dim(B)[1], dim(B)[2]), nrow = 1)
  Tn = mat.TL[1]
  Ln = mat.TL[2]
  P = dim(u)[2]
  Iden_T = Matrix(diag(Tn), sparse = T)
  Iden_T <- as.matrix(Iden_T)
  H2 = inv(t(D)%*%D)%*%t(D)
  H = D%*%H2
  ITH = Iden_T - H
  totalR = kronecker(diag(P), R)
  F_cur=obj_val(Y,D,d,B,u,v,lambda,R)
  maxIter=20
  N = dim(Y)[3] # Y_train
  for (iter in 1:maxIter)
  {
    print(iter)
    # U and d are fixed, find V 
    W = B%*%u
    Y_new <- array(data = NA, dim = c(dim(Y)[1],dim(Y)[2],dim(Y)[3]))
    for (i in 1:dim(Y)[3]) {
      Y_new[,,i] = Y[,,i]-D%*%d[,,i]
    }
    Y_new
    Q = t(W)%*%W
    inv_Q = inv(Q)
    for (i in 1:dim(Y)[3]) {
      v[,,i] = inv_Q%*%(t(W)%*%Y_new[,,i]) 
    }
    v
    # When V are fixed, 
    # min_{U} 1/n ||Y -Dd - X*U*V||^2 + lambda U' R U
    # A bit more complicated but fairly easy. 
    # Calculate sum of V_iV_i^T and XY_iV_i^T
    totalW = matrix(0, nrow = Ln*P, ncol = Ln*P)
    f = matrix(0, nrow = Ln*P, ncol = 1)
    for (n in 1:N) {
      temp = t(B)%*%ITH
      VVT = v[,,n]%*%t(v[,,n])
      XTX = temp%*%B
      totalW = totalW+ kronecker(VVT,XTX)
      f = f + matrix(temp%*%Y[,,n]%*%t(v[,,n]), nrow = P*Ln, ncol = 1)
    }
    Q = totalW/N + lambda*totalR
    f=f/N
    # pcg part
    # [sol,flag]=pcg(Q,f,1e-5,2000)
    sol <- pcg(Q, f, maxiter = 5000, tol = 1e-06)
    u = matrix(sol, nrow = Ln, ncol = P)
    # in Matlab d_cur=cell(N,1)
    # R has to set up dimension
    # For previous, d is random array with dimension (7,116,820)
    d = array(data = sol, dim = c(7,1171,1))
    for (n in 1:N) {
      d[,,n] = H2%*%Y[,,n]
      d[,,n] = d[,,n] - H2%*%B%*%u%*%v[,,n]
    }
    F_new=obj_val(Y,D,d,B,u,v,lambda,R)
    if (F_new - F_cur > 1e-5)
    {F_cur = F_new
    break}
    F_cur =  F_new
  }
  list(d=d,v=v,u=u)
}

output1 <- lr.func.stp12(one_arr, D, d, B, u, v, lambda=100, R)

output1
output1$u

matplot(output1$u, type = "l")












### test for reformat function
### m1 is 5*1 matrix
### the final output should be
###     [,1] [,2] [,3] [,4]
##[1,]    1    2    3    4
##[2,]    2    3    4    5

m1 <- matrix(data = rep(1:5), ncol = 1)
m1

reformat <- function(Y){
  h <- 2
  t <- dim(Y)[1]-h+1
  mat1 <- matrix(data = 0, nrow = h, ncol = t)
  Yn <- Y
  for (i in 1:t) {
    mat1[,i] <- Yn[i:(i+1),]
  }
  mat1
}
onep <- reformat(m1)
onep

### test for reformat function
