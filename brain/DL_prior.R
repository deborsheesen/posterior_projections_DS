rm(list=ls())
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
setwd("/home/postdoc/dsen/Desktop/projeted_posterior/brain")

model="
data{
int<lower=1> R;
int<lower=1> N;
int<lower=1> d;
int X[R,R,N];     // X is an N-array of (RxR)-matrices 
real alpha;
real lambda1;
real lambda2;
}
parameters {
matrix[R,R] Z;
matrix[R,d] U;
matrix[N,d] D;
real<lower=0> etaU[R,d];// C(0,1)
real<lower=0> sigsqZ;
real<lower=0> phiD[d];
}
model {
matrix[R,R] UDU;
matrix[d,d] U2;
vector[R] etaUSimplex;
real etaUSimplexNorm;
for (n in 1:N){
UDU = Z + U * diag_matrix(to_vector(D[n,:])) * U';
for (i in 1:R){
for (j in 1:(i-1)){
target += bernoulli_logit_lpmf(X[i,j,n] | UDU[i,j]);
}
}
}
for (s in 1:d){
target +=  normal_lpdf( to_vector(D[:,s]) | 0, phiD[s]);
target +=  cauchy_lpdf( phiD[s] | 2,1);
etaUSimplexNorm = sum(to_vector(etaU[:,s]));
etaUSimplex = to_vector(etaU[:,s])/etaUSimplexNorm;
target +=  (alpha-1)* sum(log(etaUSimplex)) ;
target +=  (R*alpha-1)*log(etaUSimplexNorm) - lambda2*  etaUSimplexNorm;
target +=  double_exponential_lpdf( to_vector(U[:,s]) | 0, etaUSimplex);
}
target +=  normal_lpdf( to_vector(Z) | 0, sqrt(sigsqZ));
target +=  gamma_lpdf( 1.0/sigsqZ | 2,1);
//orthonormality constraint  in U
U2 = U'* U - diag_matrix(rep_vector(1.0, d));
// here I used the 2-norm on U2, but it can be replaced by 1-norm
target += - lambda1* trace( U2' *U2);
}
"

load("kki42.rda")
X=tensorA[1:68,1:68,1:21]
for(i in 1:21){
  X[,,i]= t(X[,,i])
}
R = dim(X)[1]
N = dim(X)[3]
d = 10 
input.dat <- list(R=R, N=N, d=d, X=X, alpha=0.01,lambda1=100,lambda2=0.5)

#initialization
Z_ini= matrix(rnorm(R*R),R,R)
p <- apply(X, 1:2, mean)
p[p==0]<-0.01
p[p==1]<-0.99
eg_p = eigen( log(p/(1-p)))
U_ini = eg_p$vectors[,1:d]
D_ini = matrix(abs(eg_p$values[1:d]),N,d)

etaU_ini= matrix(runif(R*d),R,d)
sigsqZ_ini = runif(1)
phiD_ini = runif(d)

init.dat=list(list(
  "Z"= Z_ini,
  "U"= U_ini,
  "D"= D_ini,
  "etaU" = etaU_ini,
  "sigsqZ" = sigsqZ_ini,
  "phiD" = phiD_ini
))

fit <- stan(model_code=model, data=input.dat, init= init.dat, iter=10000, chains=1)
saveRDS(fit, file="/home/grad/sp289/Projection/fitCDL2/fitCDL2.rds")