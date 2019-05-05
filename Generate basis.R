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


