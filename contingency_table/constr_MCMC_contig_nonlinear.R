library(multinomineq)
library(dirmcmc)

model <- function(x)
  x[1]<x[6] & x[1]+x[2]<x[6]+x[7] & x[1]+x[2]+x[3]<x[6]+x[7]+x[8] &
  x[1]+x[2]+x[3]+x[4]<x[6]+x[7]+x[8]+x[9] &
  x[6]<x[11] & x[6]+x[7]<x[11]+x[12] & x[6]+x[7]+x[8]<x[11]+x[12]+x[13] &
  x[6]+x[7]+x[8]+x[9]<x[11]+x[12]+x[13]+x[14] &
  x[11]<x[16] & x[11]+x[12]<x[16]+x[17] & x[11]+x[12]+x[13]<x[16]+x[17]+x[18] &
  x[11]+x[12]+x[13]+x[14]<x[16]+x[17]+x[18]+x[19]

N = 10000
B = 2000
mcmc <- sampling_nonlinear(k=k, options=20, inside=model, M=N, burnin=B,  prior=rep(1,20))
dim(mcmc)
#
ESS = numeric(19)
for (i in 1:19) ESS[i] = (M-B)/iact(mcmc[,i])

post_samples = array(0,dim=c(dim(mcmc)[1],4,5))
for (n in 1:dim(mcmc)[1]) {
  post_samples[n,,1:4] = matrix(mcmc[n,],4,4, byrow=T)
  for (j in 1:4){
    post_samples[n,j,5] = 1-sum(post_samples[n,j,1:4])  
  }
}
dim(post_samples)
post_samples[1,,]

post_mean = apply(post_samples, c(2,3), mean)

X = matrix(k, 4,5, byrow=T) 
predicted_means_constr = array(0,dim=c(4,5))
for (i in 1:4){
  predicted_means_constr[i,] = rowSums(X)[i]*post_mean[i,]
}

error_constr = X/predicted_means_constr
mean(error_constr)