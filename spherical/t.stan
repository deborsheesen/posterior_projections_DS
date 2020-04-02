data
{
    int<lower=0> m;
    vector[m] XBar;
    vector[m] mu;
}
parameters 
{
    vector[m] theta;
}
model 
{
    for (i in 1:m)
    {
        theta[i] ~ student_t(3., mu[i], 1.);
        XBar[i] ~ normal(theta[i], sqrt(0.1));
    }
}