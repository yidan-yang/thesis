library('R.matlab')
library('matrixcalc')
setwd('C:\\Users\\yancx\\Desktop\\Thesis\\Peter fMRI data\\')
basis <- readMat("basis.mat")

# load("Scans.arr")
# Y_oneScan <- Scans.arr[,,1,]
# save(Y_oneScan, file = "Y_oneScan")

load("Y_oneScan")
Y <- Y_oneScan


# d is random array with dimension (7,116,820)
# u is random matrix with dimension (33,p) where p=2 in this case
# v is random array with dimension (p,116,820) where p=2

p <- 2
v <- array(data = runif(p*116*820), dim = c(p,116,820))
d <- array(data = runif(7*116*820), dim = c(7,116,820))
u <- matrix(data = runif(33*p), nrow = 33, ncol = p)


D <- basis$D
B <- basis$B
R <- basis$R


###################################################
#       predict_curve
###################################################
lr_lsq <- function(Y_test, D, B , U_cur){
  N <- dim(Y_test)[3]
  P <- ncol(U_cur)
  W <- B%*%U_cur
  A <- cbind(W,D) 
  inv_A = solve(t(A)%*%A) 
  # in Matlab d_test=cell(N,1), V_test=cell(N,1)
  d_test <- array(data = 0, dim = c(7,116,N))
  V_test <- array(data = 0, dim = c(2,116,N))
  for (n in 1:N) {
    sol = inv_A%*%(t(A)%*%Y_test[,,n])
    V_test[,,n] = sol[1:P,]
    D_test[,,n] = sol[-c(1:P),]
  }
  for (i in 1:dim(Y_test)[3]) {
    temp = temp + frobenius.norm(Y_test[,,i]-D%*%D_test[,,i])^2 
  }
  error = temp/N
}


###################################################
#       cross validation
###################################################

# P is a 1*19 matrix
P <- matrix(data = 2:20, nrow = 1)

# lambda is a 1*13 matrix
lambda <- matrix(data = -1:11, nrow = 1)
lambda <- exp(lambda)

error <- matrix(data = 0, nrow = length(P), ncol = length(lambda))

# cross validation
for (fold in 1:5) {
  # Split data into training and testing
  index = sample(rep(1:5, 820/5))
  train = which(index!=fold)
  test =  which(index==fold)
  Y_train = Y[,,train]
  Y_test = Y[,,test]
  ###
  # Choose the number of principal curves
  n = dim(Y_train)[3]
  T = dim(Y_train[,,1])[1]
  J = dim(Y_train[,,1])[2]
  T = dim(B)[1]
  L = dim(B)[2]
  r = dim(D)[2]
  # Initilize starting solution
  # d is random array with dimension (7,116,820)
  # u is random matrix with dimension (33,p) where p=2 in this case
  # v is random array with dimension (p,116,820) where p=2
  # here n = 820
  # here J = dim(Y_train[,,1])[2]=116
  d_cur = array(data = runif(r*J*n), dim = c(r,J,n))
  for (p in 1:length(P)) {
    U_cur = matrix(data = runif(L*P[p]), nrow = L, ncol = P[p])
    V_cur = array(data = P[p]*J*n, dim = c(P[p],J,n))
  }
  for (i in 1:length(lambda)) {
    list(d_cur=d_cur,U_cur=U_cur,V_cur=V_cur)=
      lr.func.stp12(Y_train,D,d_cur,B,U_cur,V_cur,lambda[i],R)
    list(d_test=d_test,V_test=V_test,temp=temp)=
      predict_curve(Y_test,D,B,U_cur)
    
  }
}








