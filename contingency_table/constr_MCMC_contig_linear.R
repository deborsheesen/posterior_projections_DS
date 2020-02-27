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
B = 1
samp <- sampling_multinom(k=k, options=options, prior=rep(1,20), A=-A, b=b, M=N, burnin=B)

post_samples = array(0,dim=c(4,5,dim(samp)[1]))
for (n in 1:dim(samp)[1]) {
  post_samples[,1:4,n] = matrix(samp[n,],4,4, byrow=T)
  for (j in 1:4){
    post_samples[j,5,n] = 1-sum(post_samples[j,1:4,n])  
  }
}
post_mean = apply(post_samples, c(1,2), mean)

# Bias:
predicted_means_constr = array(0,dim=c(4,5))
for (i in 1:4){predicted_means_constr[i,] = rowSums(X)[i]*post_mean[i,]}
mean(abs(X-predicted_means_constr))

par(mfrow=c(1,2),mar=c(4,0.2,2,0.5), oma=c(1, 2, 0.2, 0.2))
plot(post_samples[1,1,], type="l", xlab="", ylab="", ylim=c(0.12,0.3)) 
title(main="Traceplot for Gibbs sampler",xlab="MCMC iteration", ylab="theta_{11}")
plot(piPostProj[1,1,], type="l", xlab="", ylab="", yaxt='n') 
title(main="Samples for projections approach", xlab="Sample", ylab="", ylim=c(0.12,0.3))

CI_Gibbs = array(0, dim=c(2,4,5))
CI_proj = array(0, dim=c(2,4,5))
CI_ratio = array(0, dim=c(4,5))

for (i in 1:4){
  for (j in 1:5){
    CI_Gibbs[,i,j] = quantile(post_samples[i,j,], c(0.025,0.975))
    CI_proj[,i,j] = quantile(piPostProj[i,j,], c(0.025,0.975))
    CI_ratio[i,j] = (CI_Gibbs[2,i,j]-CI_Gibbs[1,i,j])/(CI_proj[2,i,j]-CI_proj[1,i,j])
  }
}




