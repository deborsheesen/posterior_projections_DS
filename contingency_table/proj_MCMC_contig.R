library(quadprog)

rDirichlet = function(n, alpha){
  l = length(alpha)
  x = matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  sm = rowSums(x)
  return(x/sm)
}

k = c(59,48,44,43,25, 21,14,4,46,44, 54,49,48,47,64, 58,32,30,31,41)
X = matrix(k, 4,5)  #note: **NOT** byrow=T

E = matrix(0,4,20)   # Equality restriction matrix E
for (i in 1:4) {E[i,seq(i,20,4)] = 1}

IE = matrix(0,12,20) # Inequality restriction matrix E
for (i in 1:3){
 for (j in 1:4){
   m = matrix(0,4,5)
   m[i,1:j] = 1
   m[(i+1),1:j] = -1
   ind = 4*(i-1)+j
   IE[ind,] = c(m)
 }
}
A = rbind(E, -IE) # Total matrix
bvec = c(rep(1,4),rep(0,12))


Dmat = diag(length(k))
Amat = t(A)
a = 1 


# Generate posterior samples:
N = 8000
piPost = array(0,dim=c(4,5,N))
for (i in 1:4){
  piPost[i,,] = t(rDirichlet(N, alpha = (X[i,]+a)))
}
piPostMean = apply(piPost, c(1,2), mean)

# Project posterior samples:
piPostProj = array(0,dim=c(4,5,N))
distPost = double(N)
for (n in 1:N){
  sol = solve.QP(Dmat, c(piPost[,,n]), Amat, bvec, meq = 4)$solution
  piPostProj[,,n] = matrix(sol,4,5)
  distPost[n] = norm(piPost[,,n]-piPostProj[,,n], type = "F")
}
piPostProjMean = apply(piPostProj, c(1,2), mean)

# Bias:
predicted_means = array(0,dim=c(4,5))
for (i in 1:4){predicted_means[i,] = rowSums(X)[i]*piPostMean[i,]}
mean(abs(X-predicted_means))

predicted_means_proj = array(0,dim=c(4,5))
for (i in 1:4){
  predicted_means_proj[i,] = rowSums(X)[i]*piPostProjMean[i,]
}
mean(abs(X-predicted_means_proj))









# ================================================================================================
# ================================================================================================

# OTHER STUFF

## Generate prior samples:
piPrior = array(0,dim=c(4,5,N))
for (i in 1:4){
  piPrior[i,,] = t(rDirichlet(N, alpha = rep(a,5)))
}

# Find prior mean:
piPriorMean = apply(piPrior,1:2, mean)
piPriorCum = t(apply(matrix(piPriorMean, 4, 5), 1, cumsum)) # Sanity check: Last column is 1

# Project prior samples:
piPriorProj = array(0,dim=c(4,5,N))
distPrior = double(N)
for (n in 1:N){
  temp = c(piPrior[,,n])
  sol = solve.QP(Dmat, temp, Amat, bvec, meq = 4)$solution
  piPriorProj[,,n] = matrix(sol,4,5)
  
  distPrior[n] = norm(piPrior[,,n]-piPriorProj[,,n], type = "F")
}


piPostProjLow = apply(piPostProj, 1:2, function(i)quantile(i, probs = c(.025)))
piPostProjHigh = apply(piPostProj, 1:2, function(i)quantile(i, probs = c(.975)))

piPostProjCum = t(apply(matrix(piPostProjMean,4, 5), 1, cumsum)) # Sanity check: Last column is 1


## Cumulative odds ratio

# Posterior CU
CUPost = array(0, dim=c(3,4,N))
for (n in 1:N){
  tempCU = t(apply(piPost[,,n], 1, cumsum))
  for(i in 1:3){
    CUPost[i,,n] = tempCU[i,1:4]/tempCU[(i+1),1:4]
  }
}

CUPostMean = apply(CUPost, 1:2, mean)
CUPostLow = apply(CUPost, 1:2, function(i)quantile(i, probs = c(.025)))
CUPostHigh = apply(CUPost, 1:2, function(i)quantile(i, probs = c(.975)))

# Projected CU
CUPostProj = array(0, dim=c(3,4,N))
for (n in 1:N){
  tempCU = t(apply(piPostProj[,,n], 1, cumsum))
  for(i in 1:3){
    CUPostProj[i,,n] = tempCU[i,1:4]/tempCU[(i+1),1:4]
  }
}

CUPostProjMean = apply(CUPostProj, 1:2, mean)
CUPostProjLow = apply(CUPostProj, 1:2, function(i)quantile(i, probs = c(.025)))
CUPostProjHigh = apply(CUPostProj, 1:2, function(i)quantile(i, probs = c(.975)))







# Plotting mean of $tau^2$ under a sequence of alpha values
Alpha = c(0.1, 0.5, 1, 5, 10)
distPriorMean = distPostMean = double(length(Alpha))
count = 1
for (a in Alpha){
  print(a)
  N = 10000
  
  # Prior
  piPrior = array(0,dim=c(4,5,N))
  for (i in 1:4){
    piPrior[i,,] = t(rDirichlet(N, alpha = rep(a,5)))
  }
  ## Projected probability
  piPriorProj = array(0,dim=c(4,5,N))
  distPrior = double(N)
  for (n in 1:N){
    temp = c(piPrior[,,n])
    sol = solve.QP(Dmat, temp, Amat, bvec, meq = 4)$solution
    piPriorProj[,,n] = matrix(sol,4,5)
    
    distPrior[n] = norm(piPrior[,,n]-piPriorProj[,,n], type = "F")
  }
  distPriorMean[count] = mean(distPrior)
  
  # Posterior
  piPost = array(0,dim=c(4,5,N))
  for (i in 1:4){
    piPost[i,,] = t(rDirichlet(N, alpha = (X[i,]+a)))
  }
  ## Projected probability
  piPostProj = array(0,dim=c(4,5,N))
  distPost = double(N)
  for (n in 1:N){
    temp = c(piPost[,,n])
    sol = solve.QP(Dmat, temp, Amat, bvec, meq = 4)$solution
    piPostProj[,,n] = matrix(sol,4,5)
    
    distPost[n] = norm(piPost[,,n]-piPostProj[,,n], type = "F")
  }
  distPostMean[count] = mean(distPost)
  
  count = count + 1
}

plot(Alpha, distPriorMean, type="l", col="red", lwd = 3)
lines(Alpha, distPostMean, col="green", lwd = 3)
legend("topright", c("Prior", "Posterior"), lty = c(1,1), lwd = c(3,3), col=c("red","green"))


