functions {
    // logistic function
    real logistic(real time, real A, real mu, real lambda) {
        real x;
        x = A/(1 + exp(4 * mu * (lambda - time)/A + 2));
        return x;
    }
}

data {
    int<lower=0> N; // individuals
    int<lower=1> M; // measurements
    int<lower = 0, upper = 1> sex[N];
    int<lower = 1, upper = N> ind[M];
    vector[M] time; // times
    vector[M] y; // weights
}

parameters {
  real<lower=0> sigma;

  real<lower = 0> A_0;
  real A_sex;
  real A_tilde[N];
  real<lower=0> sigma_A;
    
  real<lower = 0> mu_0;
  real mu_sex;
  real mu_tilde[N];
  real<lower=0> sigma_mu;
  
  real lambda_0;
  real lambda_sex;
  real lambda_tilde[N];
  real<lower=0> sigma_lambda;
}

transformed parameters{
  vector<lower = 0>[N] A;
  vector<lower = 0>[N] mu;
  vector[N] lambda;
  
  real A_i[N];
  real mu_i[N];
  real lambda_i[N];
  
  for(n in 1:N){
    
    A_i[n] = A_tilde[n] * sigma_A;
    mu_i[n] = mu_tilde[n] * sigma_mu;
    lambda_i[n] = lambda_tilde[n] * sigma_lambda;
    
    A[n]      = 10 * (A_0 + sex[n] *      A_sex +      A_i[n]);
    mu[n]     =      mu_0 + sex[n] *     mu_sex +     mu_i[n];
    lambda[n] =  lambda_0 + sex[n] * lambda_sex + lambda_i[n];
  }
}

model {
  vector[M] x;
  
  for (i in 1:(M)){
      x[i] = logistic(time[i], A[ind[i]], mu[ind[i]], lambda[ind[i]]);
  }
  y ~ normal(x, sigma);
  
  A_0 ~ normal(0.0, 5.0);
  mu_0 ~ normal(0.0, 1.0);
  lambda_0 ~ normal(0.0, 5.0);
  
  A_sex ~ normal(0.0, 1.0);
  mu_sex ~ normal(0.0, 1.0);
  lambda_sex ~ normal(0.0, 1.0);
  
  sigma ~ cauchy(0, 2.5);
  
  sigma_A ~ cauchy(0, 1);
  sigma_mu ~ cauchy(0, 1);
  sigma_lambda ~ cauchy(0, 1);

  A_tilde ~ normal(0, 1);
  mu_tilde ~ normal(0, 1);
  lambda_tilde ~ normal(0, 1);
}
