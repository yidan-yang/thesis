# output.v
setwd('C:\\Users\\yancx\\Desktop\\Thesis\\Peter fMRI data\\')
load("output.v")

class(output.v)
dim(output.v)


extract.116.v <- function(Y){
  out.arr <- array(data = 0, dim = c(100,1,2,116))
  for (j in 1:116) {
    for (i in 1:2) {
      out.arr[,,i,j] <- Y[i,(1171*(j-1)+1),]
    }
  }
  out.arr
}
out.array <- extract.116.v(output.v)


out.df <- data.frame(data=out.array)
person <- paste(rep(1:25, each=4), sep = "")
out.df <- cbind(person, out.df)


colnames(out.df) <- c("person", 1:232)


adjsq.func.v <- function(data){
  out.var <- matrix(data = 0, nrow = 232)
  for (i in 2:233) {
    out.var[(i-1),] <- summary(lm(out.df[,i] ~ person, data = out.df))$adj.r.squared
  }
  out.var <- matrix(data = out.var, nrow = 2, ncol = 116)
  out.var
}
adjsq.out.v <- adjsq.func.v(out.df)

adjsq.out.v [adjsq.out.v  < 0] <- 0
adjsq.out.v 


p.value.func <- function(data){ 
  out.pvalue <- matrix(data = 0, nrow = 232)
  for (i in 2:233) {
    out.pvalue[(i-1),] <- anova(lm(out.df[,i] ~ person, data = out.df))$'Pr(>F)'[1]
  }
  out.pvalue <- matrix(data = out.pvalue, nrow = 2, ncol = 116)
  out.pvalue
}
p.value.v <- p.value.func(out.df)
p.value.v
