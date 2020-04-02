data
{
    int<lower=1> m;
    int<lower=1> N;
    int<lower=1> p;
    int X[N,m,m];     // X is an N-array of (mxm)-matrices 
}

parameters 
{
    matrix[m,m] Z;
    matrix[m,p] U;
    matrix[N,p] D;  
    real<lower=0> sigsqZ;
    real<lower=0> sigsqD[p];
}

model 
{
    matrix[m,m] UDU;
    for (n in 1:N)
    {
        UDU = Z + U * diag_matrix(to_vector(D[n,:])) * U';
        for (i in 1:m)
        {
            for (j in 1:(i-1))
            {
                target += bernoulli_logit_lpmf(X[n,i,j] | UDU[i,j]);
            }
        }
    }
    for (s in 1:p)
    {
        target +=  normal_lpdf( to_vector(D[:,s]) | 0, sqrt(sigsqD[s]));
        target +=  gamma_lpdf( 1.0/sigsqD[s] | 2,1);
        target +=  normal_lpdf( to_vector(U[:,s]) | 0, 1);
    }
    target +=  normal_lpdf( to_vector(Z) | 0, sqrt(sigsqZ));
    target +=  gamma_lpdf( 1.0/sigsqZ | 2,1);
}