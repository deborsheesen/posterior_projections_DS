library(truncnorm)
library(latex2exp)
set.seed(3)

n = 50
theta0 = -0.5
x.bar <- mean(rnorm(n,theta0,1))
#x.bar = 0.5
#x.bar = 0, .75
#x0 <- x.bar
sigma.n.sq <- 1/(n+0.001)
theta.n <- sigma.n.sq * n * x.bar
alpha <- -(theta.n/ sqrt(sigma.n.sq))


W1 <- pnorm(alpha)
W2 <- 1 - W1
C1 <- dnorm(0, x.bar, sqrt(1/n))
C2 <- (C1 * dtruncnorm(0, 0, Inf, 0, sqrt(1000)))/dtruncnorm(0, 0, Inf, theta.n, sqrt(sigma.n.sq))

w1 <- W1/C1
w2 <- W2/C2

w1 <- w1/(sum(w1+w2))
w2 <- 1 - w1
print(w1)

Z <- 1-pnorm(alpha)
unconstrained.mean <- theta.n
trunc.mean <- theta.n + sqrt(sigma.n.sq)*dnorm(alpha)/Z
proj.mean <- Z*trunc.mean
projnorm <- function(x){
  if (x == 0) return(pnorm(alpha))
  if (x>0) return(Z*dtruncnorm(x, 0, Inf, theta.n, sqrt(sigma.n.sq)))
  else return(0)
}

prior <- function(x, w1, w2){
  if (x == 0) return(w1)
  if (x>0) return(w2*dtruncnorm(x, 0, Inf, 0, sqrt(1000)))
  else return(0)
}
a <- -0.2
b <- 1

x <- seq(a,b,.01)
x.neg <- seq(a,-0.01,0.01)
x.pos <- seq(0.01, b, 0.01)

# jpeg("eg1_4.jpeg", width = 5, height = 5, units = 'in', res = 300)


plot(x.pos, sapply(x.pos, function(ind) projnorm(ind)), col="red", type="l", 
     xlab=TeX("$\\theta$"), lwd = 3, xlim=c(-0.05,1), ylim = c(0,3), ylab="Density",
     cex.lab=1.5, cex.axis=1.3, cex.main=1.5, cex.sub=1.5)
lines(x.pos, sapply(x.pos, function(ind) prior(ind, w1, w2)), col="green", 
      lwd=3, lty=2)
grid()
segments(0,0,0,projnorm(0),col="red",lwd=3)
segments(0,0,0,prior(0,w1,w2),col="green",lwd=3, lty = 2)
abline(h=0, lwd=2,col=rgb(0,0,0,alpha=0.5))

legend(0.55,3,c("Prior", "Posterior"), bty="n",
       col=c("green","red"),lty=c(1,2), lwd=3,cex=1.4)

#dev.off()