functions {
    // gompertz function
    real gompertz(real time, real A, real mu, real lambda) {
        real x;
        x = A * exp(-exp(mu * exp(1) * (lambda - time)/A + 1));
        return x;
    }
}

data {
    int<lower=0> N; // individuals
    int<lower=1> K; // traits
    int<lower = 1, upper = N> ind[K*N];
    vector[K*N] time; // times
    vector[K*N] y; // weights
}

parameters {
  real<lower = 0> A_0;
  real<lower = 0> mu_0;
  real lambda_0;
  real sigma;
  
  real sigma_A;
  real sigma_mu;
  real sigma_lambda;
  
  vector[N] A_i;
  vector[N] mu_i;
  vector[N] lambda_i;
}

model {
  vector[K*N] x;
  
  for (i in 1:(N*K)){
      x[i] = gompertz(time[i], 10 * A_0, mu_0, lambda_0);
  }
  y ~ normal(x, sigma);
  
  A_0 ~ normal(0.0, 5.0);
  mu_0 ~ normal(0.0, 1.0);
  lambda_0 ~ normal(0.0, 5.0);
  
  sigma ~ cauchy(0, 2.5);
}
