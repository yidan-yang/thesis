library('R.matlab')
library('matrixcalc')
setwd('C:\\Users\\yancx\\Desktop\\Thesis\\Peter fMRI data\\')
basis <- readMat("basis.mat")
load("Scans.arr")



Y <- Scans.arr[,,1,]

# d is random array with dimension (7,116,820)
# u is random matrix with dimension (33,p) where p=2 in this case
# v is random array with dimension (p,116,820) where p=2

p <- 2
v <- array(data = runif(p*116*820), dim = c(p,116,820))

set.seed(10)
d <- array(data = runif(7*116*820), dim = c(7,116,820))
u <- matrix(data = runif(33*p), nrow = 33, ncol = p)


D <- basis$D
B <- basis$B
R <- basis$R



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
  ########################################
  # May 6 
  # the following for loop is for rest of function
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
    #}
    ################################# iter for loop should not end here
    
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
    d = array(data = sol, dim = c(7,116,820))
    for (n in 1:N) {
      d[,,n] = H2%*%Y[,,n]
      d[,,n] = d[,,n] - H2%*%B%*%u%*%v[,,n]
    }
    F_new=obj_val(Y,D,d,B,u,v,lambda,R)
    print(F_new)
    if (F_new - F_cur > 1e-5)
    {F_cur = F_new
    break}
    F_cur =  F_new
  }
  list(d=d,v=v,u=u)
}

output1 <- lr.func.stp12(Y, D, d, B, u, v, lambda=100, R)

matplot(B%*%output1$u, type = "l", ylab = "B*u")

# d is random array with dimension (7,116,820)
# out.d is matrix with dimension (7,116*820)
out.d <- matrix(data = output1$d, ncol = 116*820)
image(out.d)

matplot(B,type = "l")
matplot(D,type = "l")

