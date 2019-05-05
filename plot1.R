library('R.matlab')
library('matrixcalc')
setwd('C:\\Users\\yancx\\Desktop\\Thesis\\Peter fMRI data\\')
basis <- readMat("basis.mat")
load("Scans.arr")
scan <- Scans.arr

dim(scan)


Y <- scan[,,1,]


D <- basis$D
B <- basis$B
R <- basis$R

dim(D)

#### plot B, dim(B)=1200*33
matplot(B, type = "l")

## plot B, h=30
## Bn is an array w/ dim 30*33*1171

reformatB <- function(B){
  h <- 30
  k <- dim(B)[2]
  t <- dim(B)[1]-h+1
  aa1 <- array(data = 0,dim = c(h,k,t))
  for (i in 1:t) {
    aa1[,,i] <- B[i:(i+29),]
  }
  aa1
}
Bn <- reformatB(B)

windows()
par(mfrow=c(3,1))
matplot(Bn[,,1], type = "l")
matplot(Bn[,,2], type = "l")
matplot(Bn[,,3], type = "l")



#### plot D, dim(D) = 1200*7
matplot(D, type = "l")

## plot D, h=30
## Dn is an array w/ dim 30*7*1171

reformatD <- function(D){
  h <- 30
  k <- dim(D)[2]
  t <- dim(D)[1]-h+1
  aa1 <- array(data = 0,dim = c(h,k,t))
  for (i in 1:t) {
    aa1[,,i] <- D[i:(i+29),]
  }
  aa1
}
Dn <- reformatD(D)

windows()
par(mfrow=c(3,1))
matplot(Dn[,,1], type = "l")
matplot(Dn[,,2], type = "l")
matplot(Dn[,,3], type = "l")


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
  aa1
}
onep <- reformat(Y)

dim(onep)


p <- 2
v <- array(data = runif(p*116*820), dim = c(p,116,820))
d <- array(data = runif(7*116*820), dim = c(7,116,820))
u <- matrix(data = runif(33*p), nrow = 33, ncol = p)



Ysub1 <- Y[,,1]
Yss <- array(data = Ysub1, dim = c(1200,116,1))



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
obj_val(Yss,D,d,B,u,v,lambda = 100,R)
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
outputsub1 <- lr.func.stp12(Yss, D, d, B, u, v, lambda=100, R)

matplot(outputsub1$u, type = "l")





