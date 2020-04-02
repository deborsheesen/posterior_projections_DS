library(Directional)
library(MASS)
library(mvtnorm)

set.seed(3)


## Generate data
n = 100
#theta = rep(1/sqrt(3), 3)
theta = c(4,2,1)
theta = theta/sqrt(sum(theta^2))
X = mvrnorm(n, mu=theta, Sigma = 10*diag(3))
XBar = colMeans(X)

#################################################
## Posterior : VMF
#################################################
mu = rep(1/sqrt(3), 3)
#mu = theta
mu = mu/sqrt(sum(mu^2))
phi = 1
temp = 0.1*n*XBar + phi*mu
phin = sqrt(sum(temp^2))
mun = temp/ phin
x = rvmf(1000, mun, phin)
dat = data.frame("x" = x[,1], "y" = x[,2], "z" =x[,3])

# jpeg("eg2FisherPost.jpeg", width = 5, height = 5, units = 'in', res = 300)
pairs(dat, main="von Mises-Fisher distribution")
#dev.off()
dat.mean <- colMeans(dat)
mean((dat.mean - theta)^2)
apply(dat, 2, function(x) quantile(x,probs = c(0.025, 0.975)))


#################################################
## Posterior : Normal
#################################################
sigman <- 1/(phi + 0.1*n)
mun <- sigman*temp
y <- mvrnorm(1000, mu=mun, Sigma = sigman*diag(3))
yProj <- y
for (i in 1:1000){
  yProj[i,] <- y[i,]/sqrt(sum(y[i,]^2))
}
dat2 <- data.frame("x" = yProj[,1], "y" = yProj[,2], "z" =yProj[,3])


# jpeg("eg2NormPost.jpeg", width = 5, height = 5, units = 'in', res = 300)
pairs(dat2, main = "Projected Spherical Normal distribution")
#dev.off()
dat2.mean = colMeans(dat2)
mean((dat2.mean - theta)^2)
apply(dat2, 2, function(x) quantile(x,probs = c(0.025, 0.975)))
#################################################
## Posterior : t
#################################################
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')

model="
data{
real XBar[3];
real mu[3];
}
parameters {
real theta[3];
}
model {
//Prior
theta ~ student_t(3, mu, sqrt(0.1));
//Data
XBar ~ normal(theta, sqrt(0.1));
}
"
input.dat = list(XBar = XBar, mu = mu)
fit = stan(model_code=model, data=input.dat, iter=2000, chains=1)
postFit = extract(fit,permuted=TRUE)
z = postFit$theta

zProj = z
for (i in 1:1000){
  zProj[i,] = z[i,]/sqrt(sum(z[i,]^2))
}
dat3 = data.frame("x" = zProj[,1], "y" = zProj[,2], "z" =zProj[,3])

# jpeg("eg2tPost.jpeg", width = 5, height = 5, units = 'in', res = 300)
pairs(dat3, main = "Projected Sperical t distribution")
# dev.off()