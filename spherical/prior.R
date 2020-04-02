library(Directional)
library(MASS)
library(mvtnorm)

set.seed(3)

m <- rep(1, 3)
m <- c(2,1,0)
m <- m/sqrt(sum(m^2))
x <- rvmf(1000, m, 10)
dat <- data.frame("x" = x[,1], "y" = x[,2], "z" =x[,3])
pairs(dat, main="von Mises-Fisher distribution")

y <- mvrnorm(1000, mu=m, Sigma = 0.1*diag(3))
yProj <- y
for (i in 1:1000){
  yProj[i,] <- y[i,]/sqrt(sum(y[i,]^2))
}
dat2 <- data.frame("x" = yProj[,1], "y" = yProj[,2], "z" =yProj[,3])
pairs(dat2, main = "Projected Spherical Normal distribution")

z <- rmvt(1000, delta=m, sigma = 0.1*diag(3), df=3)
zProj <- z
for (i in 1:1000){
  zProj[i,] <- z[i,]/sqrt(sum(z[i,]^2))
}
dat3 <- data.frame("x" = zProj[,1], "y" = zProj[,2], "z" =zProj[,3])
pairs(dat3, main = "Projected Sperical t distribution")

