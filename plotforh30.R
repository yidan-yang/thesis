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
v <- array(data = runif(p*116), dim = c(p,116,1))
d <- array(data = runif(7*116), dim = c(7,116,1))
u <- matrix(data = runif(33*p), nrow = 33, ncol = p)

# For patient 1
Y1 <- Y[,,1]

# For patient1 subset(h=30) 1-30
# matrix w/ dim 30*116

Y11 <- Y1[1:30,]
Y11 <- array(data = Y11, dim = c(30,116,1))




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
obj_val(Y11,D,d,B,u,v,lambda = 100,R)
###################################################

# step1: fix U and d, solve v(v is unknown here)
# step2: fix v(at solution of previous step), u and d are unknown here

library("Matrix") # sparse matrix
library("pcg")

lr.func.stp12 <- function(Y, D, d, B, u, v, lambda, R){
  # [T,L]=size(X)
  mat.TL <- matrix(c(dim(B)[1], dim(B)[2]), nrow = 1)
  Tn = mat.TL[1]
  Ln = mat.TL[2]
  P = dim(u)[2]
  Iden_T = Matrix(diag(Tn), sparse = T)
  Iden_T <- as.matrix(Iden_T)
  H2 = solve(t(D)%*%D)%*%t(D)
  H = D%*%H2
  ITH = Iden_T - H
  totalR = kronecker(diag(P), R)
  F_cur=obj_val(Y,D,d,B,u,v,lambda,R)
  maxIter=20
    N = dim(Y)[3] # Y_train
   for (iter in 1:maxIter) {
    # U and d are fixed, find V 
    W = B%*%u
    Y_new <- array(data = NA, dim = c(dim(Y)[1],dim(Y)[2],dim(Y)[3]))
    for (i in 1:dim(Y)[3]) {
      Y_new[,,i] = Y[,,i]-D%*%d[,,i]
    }
    Y_new
    Q = t(W)%*%W
    inv_Q = solve(Q)
    for (i in 1:dim(Y)[3]) {
      v[,,i] = inv_Q%*%(t(W)%*%Y_new[,,i]) 
    }
    v
  }
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
  sol <- pcg(Q, f, maxiter = 2000, tol = 1e-05)
  u = matrix(sol, nrow = Ln, ncol = P)
  # in Matlab d_cur=cell(N,1)
  # R has to set up dimension
  # For previous, d is random array with dimension (7,116,820)
  d = array(data = sol, dim = c(7,116,820))
  for (n in 1:N) {
    d[,,n] = H2%*%Y[,,n]
    d[,,n] = d[,,n] - H2%*%B%*%u%*%v[,,n]
  }
  F_new=obj_val(Y,D,d,B,u,v,lambda,R)
  if (F_new - F_cur > 1e-5){F_cur = F_new}
  F_cur =  F_new
  list(d=d,v=v,u=u)
}

### h=30, window1
outputsub1 <- lr.func.stp12(Y11, D, d, B, u, v, lambda=100, R)
matplot(outputsub1$u, type = "l")

### h=30, window2
Y12 <- Y1[2:31,]
Y12 <- array(data = Y12, dim = c(30,116,1))

outputsub2 <- lr.func.stp12(Y12, D, d, B, u, v, lambda=100, R)
matplot(outputsub2$u, type = "l")

### h=30, window3
Y13 <- Y1[3:32,]
Y13 <- array(data = Y13, dim = c(30,116,1))

outputsub3 <- lr.func.stp12(Y13, D, d, B, u, v, lambda=100, R)
matplot(outputsub3$u, type = "l")

### h=30, window4
Y14 <- Y1[4:33,]
Y14 <- array(data = Y14, dim = c(30,116,1))

outputsub4 <- lr.func.stp12(Y14, D, d, B, u, v, lambda=100, R)
matplot(outputsub4$u, type = "l")

windows()
par(mfrow=c(4,1))
matplot(outputsub1$u, type = "l")
matplot(outputsub2$u, type = "l")
matplot(outputsub3$u, type = "l")
matplot(outputsub4$u, type = "l")
