data
{
    int<lower=1> m;
    int<lower=1> N;
    int<lower=1> p;
    int X[N,m,m];     // X is an N-array of (mxm)-matrices 
    real alpha;
    real<lower=0> tau;
}

parameters 
{
    matrix[m,m] Z;
    matrix[m,p] U;
    matrix[N,p] D;  // this is Lambda in the paper
    real<lower=0> etaU[m,p];// C(0,1)
    real<lower=0> sigsqZ;
    real<lower=0> sigsqD[p];
}

model 
{
    matrix[m,m] UDU;
    real lmbda;
    matrix[p,p] U2;
    vector[m] etaUSimplex;
    real etaUSimplexNorm;
    
    lmbda = 1/2;
    
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
        target +=  normal_lpdf(to_vector(D[:,s]) | 0, sqrt(sigsqD[s]));
        target +=  gamma_lpdf(1.0/sigsqD[s] | 2,1);
        etaUSimplexNorm = sum(to_vector(etaU[:,s]));
        etaUSimplex = to_vector(etaU[:,s])/etaUSimplexNorm;
        target +=  (alpha-1)* sum(log(etaUSimplex));  # Dirichlep distribution
        target +=  (m*alpha-1)*log(etaUSimplexNorm) - lmbda*etaUSimplexNorm;
        target +=  double_exponential_lpdf(to_vector(U[:,s]) | 0, etaUSimplex);
    }
    target +=  normal_lpdf(to_vector(Z) | 0, sqrt(sigsqZ));
    target +=  gamma_lpdf(1.0/sigsqZ | 2,1);
    
    //orthonormality constraint  in U
    U2 = U'*U - diag_matrix(rep_vector(1.0, p));
    target += - 1/tau * trace(U2'*U2);
}

