data
{
    int<lower=1> R;
    int<lower=1> N;
    int<lower=1> d;
    int X[N,R,R];     // X is an N-array of (RxR)-matrices 
    real alpha;
    real lmbda;
}

parameters 
{
    matrix[R,R] Z;
    matrix[R,d] U;
    matrix[N,d] D;
    real<lower=0> etaU[R,d];// C(0,1)
    real<lower=0> sigsqZ;
    real<lower=0> sigsqD[d];
}

model 
{
    matrix[R,R] UDU;
    vector[R] etaUSimplex;
    real etaUSimplexNorm;
    
    for (n in 1:N)
    {
        UDU = Z + U * diag_matrix(to_vector(D[n,:])) * U';
        for (i in 1:R)
        {
            for (j in 1:(i-1))
            {
                target += bernoulli_logit_lpmf(X[n,i,j] | UDU[i,j]);
            }
        }
    }

    for (s in 1:d)
    {
        target +=  normal_lpdf( to_vector(D[:,s]) | 0, sqrt(sigsqD[s]));
        target +=  gamma_lpdf( 1.0/sigsqD[s] | 2,1);
        etaUSimplexNorm = sum(to_vector(etaU[:,s]));
        etaUSimplex = to_vector(etaU[:,s])/etaUSimplexNorm;
        target +=  (alpha-1)* sum(log(etaUSimplex)) ;
        target +=  (R*alpha-1)*log(etaUSimplexNorm) - lmbda*  etaUSimplexNorm;
        target +=  double_exponential_lpdf( to_vector(U[:,s]) | 0, etaUSimplex);
    }
    target +=  normal_lpdf( to_vector(Z) | 0, sqrt(sigsqZ));
    target +=  gamma_lpdf( 1.0/sigsqZ | 2,1);
}