a.arr <- array(data = rep(1:144), dim = c(6,3,4,2))
a.arr



reformat4scan <- function(Y){
  h <- 3
  k <- dim(Y)[2]
  t <- dim(Y)[1]-h+1
  sn <- dim(Y)[3]
  subn <- dim(Y)[4]
  fin.out.arr <- array(data = 0, dim = c(h,t*k,sn,subn))
  for (index in 1:sn) {
    # aa1 is for one patient MA time series
    # aa1 will eventually convert into matrix then store in an array aa2
    # finally it will convert into matrix w/ dim(h,t*k)
    aa1 <- array(data = 0,dim = c(h,k,t))
    # aa2 is for all 820 patients
    # also the final output array w/ dim 30*(1171*116)*820
    aa2 <- array(data = 0,dim = c(h,t*k,subn))
    # extract one scan
    Y1s <- Y[,,index,]
      for (ind in 1:subn) {
        for (i in 1:t) {
          aa1[,,i] <- Y1s[i:(i+2),,ind]
          aa2[,,ind] <- matrix(aa1,ncol = t*k)
        }
        
      }
    fin.out.arr[,,index,] <- aa2
  }
  fin.out.arr <- array(data = fin.out.arr, dim = c(h,t*k,sn*subn))
  fin.out.arr
}
reYn <- reformat4scan(a.arr)
reYn

dim(reYn)
a.arr


#############################################################



reformat2 <- function(Y){
  h <- 3
  k <- dim(Y)[2]
  t <- dim(Y)[1]-h+1
  # aa1 is for one patient MA time series
  # aa1 will eventually convert into matrix then store in an array aa2
  # finally it will convert into matrix w/ dim(h,t*k)
  aa1 <- array(data = 0,dim = c(h,k,t))
  # aa2 is for all 820 patients
  # also the final output array w/ dim 30*(1171*116)*820
  aa2 <- array(data = 0,dim = c(h,t*k,dim(Y)[3]))
  for (ind in 1:dim(Y)[3]) {
    for (i in 1:t) {
      aa1[,,i] <- Y[i:(i+2),,ind]
      aa2[,,ind] <- matrix(aa1,ncol = t*k)
    }
  }
  aa2
}
reYn <- reformat2(a.test)
reYn
a.arr


