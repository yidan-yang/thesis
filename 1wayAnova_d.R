# output.d

setwd('C:\\Users\\yancx\\Desktop\\Thesis\\Peter fMRI data\\')
load("output.d")

class(output.d)
dim(output.d)

Y1 <- output.d[1,1,]
Y1 <- matrix(data = Y1 , ncol =1)
Y2 <- output.d[2,1,]
Y2 <- matrix(data = Y2 , ncol =1)
Y3 <- output.d[3,1,]
Y3 <- matrix(data = Y3 , ncol =1)
Y4 <- output.d[4,1,]
Y4 <- matrix(data = Y4 , ncol =1)
Y5 <- output.d[5,1,]
Y5 <- matrix(data = Y5 , ncol =1)
Y6 <- output.d[6,1,]
Y6 <- matrix(data = Y6 , ncol =1)
Y7 <- output.d[7,1,]
Y7 <- matrix(data = Y7 , ncol =1)


extract.116.test <- function(Y){
  out.arr <- array(data = 0, dim = c(100,1,7,116))
  for (j in 1:116) {
    for (i in 1:7) {
      out.arr[,,i,j] <- Y[i,(1171*(j-1)+1),]
    }
  }
  out.arr
}
out.array <- extract.116.test(output.d)


out.df <- data.frame(data=out.array)
person <- paste(rep(1:25, each=4), sep = "")
out.df <- cbind(person, out.df)


colnames(out.df) <- c("person", 1:812)


adjsq.func.d <- function(data){
  out.var <- matrix(data = 0, nrow = 812)
  for (i in 2:813) {
  out.var[(i-1),] <- summary(lm(out.df[,i] ~ person, data = out.df))$adj.r.squared
  }
  out.var <- matrix(data = out.var, nrow = 7, ncol = 116)
  out.var
}
adjsq.out.d <- adjsq.func.d(out.df)

adjsq.out.d [adjsq.out.d  < 0] <- 0
adjsq.out.d 


p.value.func <- function(data){ 
  out.pvalue <- matrix(data = 0, nrow = 812)
  for (i in 2:813) {
    out.pvalue[(i-1),] <- anova(lm(out.df[,i] ~ person, data = out.df))$'Pr(>F)'[1]
  }
  out.pvalue <- matrix(data = out.pvalue, nrow = 7, ncol = 116)
  out.pvalue
}
p.value.d <- p.value.func(out.df)
p.value.d



# fit <- lm(out.df[,6] ~ person, data = out.df)
# extract p-value from model
# anova(fit)$'Pr(>F)'[1]
# summary(fit)$adj.r.squared
# anova(lm(out.df[,2] ~ person, data = out.df))["person", "Mean Sq"]

# head(out.df)[1:6]



#######################################################
#           test
#######################################################
t1 <- array(data = rep(1:72), dim = c(2,3,4,3))
t11 <- array(data = t1, dim = c(2,12,3))
t11[1,1,]


extract.116.test <- function(Y){
  out.arr <- array(data = 0, dim = c(3,1,2,4))
  for (j in 1:4) {
      for (i in 1:2) {
    out.arr[,,i,j] <- Y[i,(3*(j-1)+1),]
      }
  }
  out.arr
}
extract.116.test(t11)



a1 <- array(data = rep(1:24),dim = c(3,1,2,4))
a11 <- array(data = a1, dim = c(3,1,8))
a1.df <- data.frame(data=a1)









