library(multinomineq)
library(dirmcmc)

k = c(59,48,44,43,25, 21,14,4,46,44, 54,49,48,47,64, 58,32,30,31,41)
X = matrix(k, 4,5)  #note: **NOT** byrow=T

A = matrix(c(-1, 0, 0, 0,  1, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
             -1,-1, 0, 0,  1, 1, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
             -1,-1,-1, 0,  1, 1, 1, 0,  0, 0, 0, 0,  0, 0, 0, 0,
             -1,-1,-1,-1,  1, 1, 1, 1,  0, 0, 0, 0,  0, 0, 0, 0,
              0, 0, 0, 0, -1, 0, 0, 0,  1, 0, 0, 0,  0, 0, 0, 0,
              0, 0, 0, 0, -1,-1, 0, 0,  1, 1, 0, 0,  0, 0, 0, 0,
              0, 0, 0, 0, -1,-1,-1, 0,  1, 1, 1, 0,  0, 0, 0, 0,
              0, 0, 0, 0, -1,-1,-1,-1,  1, 1, 1, 1,  0, 0, 0, 0,
              0, 0, 0, 0,  0, 0, 0, 0, -1, 0, 0, 0,  1, 0, 0, 0,
              0, 0, 0, 0,  0, 0, 0, 0, -1,-1, 0, 0,  1, 1, 0, 0, 
              0, 0, 0, 0,  0, 0, 0, 0, -1,-1,-1, 0,  1, 1, 1, 0,
              0, 0, 0, 0,  0, 0, 0, 0, -1,-1,-1,-1,  1, 1, 1, 1),
              ncol=16, byrow=TRUE)
  
b = c(0,0,0,0,0,0,0,0,0,0,0,0)

options = c(5,5,5,5)

N = 10000
B = 2000
samp <- sampling_multinom(k=k, options=options, prior=rep(1,20), A=-A, b=b, M=N, burnin=B)

post_samples = array(0,dim=c(dim(samp)[1],4,5))
for (n in 1:dim(samp)[1]) {
  post_samples[n,,1:4] = matrix(samp[n,],4,4, byrow=T)
  for (j in 1:4){
    post_samples[n,j,5] = 1-sum(post_samples[n,j,1:4])  
  }
}
post_mean = apply(post_samples, c(2,3), mean)

# Bias:
predicted_means_constr = array(0,dim=c(4,5))
for (i in 1:4){predicted_means_constr[i,] = rowSums(X)[i]*post_mean[i,]}
mean(abs(X-predicted_means_constr))

