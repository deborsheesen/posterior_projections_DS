
install.packages("multinomineq")
library(multinomineq)


################################# BAYES FACTOR

?bf_nonlinear

##### 2x2x2 continceny table (Klugkist & Hojtink, 2007)
#
# (defendant's race) x (victim's race) x (death penalty)
# indexing: 0 = white/white/yes  ; 1 = black/black/no
# probabilities: (p000,p001,  p010,p011,  p100,p101,  p110,p111)
# Model2:
# p000*p101 < p100*p001  &   p010*p111 < p110*p011

# observed frequencies:
k <- c(19,132,   0,9,   11,52,   6,97)

model <- function(x)
  x[1]*x[6] < x[5]*x[2]  &  x[3]*(1-sum(x)) < x[7]*x[4]
# NOTE: "1-sum(x)"  must be used instead of "x[8]"!

# compute Bayes factor (Klugkist 2007: bf_0u=1.62)
bf_nonlinear(k = k, options = 8, inside = model, M = 50000)





################################# PARAMETER ESTIMATION / MCMC SAMPLING

# parameter estimation: only for convex constraints (does not work for the example above!)
?sampling_nonlinear

# two binomial success probabilities: x = c(x1, x2)
# restriction to a circle:
model <- function(x)
  (x[1]-.50)^2 + (x[2]-.50)^2 <= .15

# draw prior samples
mcmc <- sampling_nonlinear(k = 0, options = c(2,2),
                           inside = model, M = 1000)
head(mcmc)
plot(c(mcmc[,1]), c(mcmc[,2]), xlim=0:1, ylim=0:1)



################################# PARAMETER ESTIMATION / REJECTION SAMPLING

# see first example above: contengency table
k <- c(19,132,   0,9,   11,52,   6,97)

model <- function(x)
  x[1]*x[6] < x[5]*x[2]  &  x[3]*(1-sum(x)) < x[7]*x[4]

# 1. sample multinomial probabilities from the unconstrained posterior
#    for this purpose: define a dummy constraint: 0*p1 + 0*p2 + ... + 0*p7 < 1
A = matrix(rep(0,8-1), nrow = 1)
b = 1
mcmc_unconstrained <- sampling_multinom(k, options = 8, A=A, b=b, M = 1e5)
head(mcmc_p)

# 2. reject all samples that violate constraints
constraint_satisfied <- apply(mcmc_p, 1, model)
mcmc_constrained <- mcmc_unconstrained[constraint_satisfied,]
nrow(mcmc_constrained)  # number of retained samples

# 3. get posterior summary
summary(mcmc_constrained)

# compare to unconstrained model:
summary(mcmc_unconstrained)
