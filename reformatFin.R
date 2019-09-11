# setwd('C:\\Users\\yancx\\Desktop\\Thesis\\Peter fMRI data\\')
setwd('/home/yy2513/fMRIdata/')

load("Scans.arr")

Y_25 <- Scans.arr[,,,1:25]


reformat4scan <- function(Y){
  h <- 30
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
        aa1[,,i] <- Y1s[i:(i+29),,ind]
        aa2[,,ind] <- matrix(aa1,ncol = t*k)
      }
    }
    fin.out.arr[,,index,] <- aa2
  }
  fin.out.arr <- array(data = fin.out.arr, dim = c(h,t*k,sn*subn))
  fin.out.arr
}
Scan30_25.arr <- reformat4scan(Y_25)

save(Scan30_25.arr, file = "Scan30_25.arr")

dim(Scan30_25.arr)

