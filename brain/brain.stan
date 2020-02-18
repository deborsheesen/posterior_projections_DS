data 
{
    int<lower=0> N;
    int<lower=0> p;
    int<lower=0> m;
    int<lower=0,upper=1> X[N,m,m];
}

parameters 
{
    matrix[m,m] Z;
    vector[p] u[m];
    vector[p] Lambda[N];
    real<lower=0> sigma_Z2;
    real<lower=0> sigma2[p];
    real<lower=0> sigma_U2;
}

model 
{
    
    for (i in 1:m) 
    {
        Z[:,i] ~ normal(0, sigma_Z2);
    }
    for (n in 1:N) 
    {
        for (k in 1:p) 
        {
            Lambda[n][k] ~ normal(0, sigma2[k]);
        }
    }
    sigma_Z2 ~ inv_gamma(2,1);
    for (k in 1:p) 
    {
        sigma2[k] ~ inv_gamma(2,1); 
    }
    for (i in 1:m) 
    {
        for (k in 1:p) 
        {
            u[i][k] ~ normal(0, sigma_U2);
        }
    }
    
    for (n in 1:N)
    {
        for (i in 1:m) 
        {
            for (j in 1:m) 
            {
                X[n,i,j] ~ bernoulli_logit(Z[i,j] + u[i]'*diag_matrix(Lambda[n])*u[j]);
            }
        }
    }
}

