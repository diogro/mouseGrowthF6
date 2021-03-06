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
    vector<lower = 0, upper = 1>[N] sex;
    int<lower = 1, upper = N> ind[M];
    vector[M] y; // weights
    int<lower = 0, upper = 1> run_estimation;
    vector[M] time; // times
}

parameters {
  real<lower=0> sigma;

  real<lower = 0> A_0;
  real A_sex;
  vector[N] A_tilde;
  real<lower=0> sigma_A;
    
  real<lower = 0> mu_0;
  real mu_sex;
  vector[N] mu_tilde;
  real<lower=0> sigma_mu;
  
  real lambda_0;
  real lambda_sex;
  vector[N] lambda_tilde;
  real<lower=0> sigma_lambda;
}

transformed parameters{
  vector[N] A;
  vector[N] mu;
  vector[N] lambda;
  
  vector[N] A_i;
  vector[N] mu_i;
  vector[N] lambda_i;
  
  A_i = A_tilde * sigma_A;
  mu_i = mu_tilde * sigma_mu;
  lambda_i = lambda_tilde * sigma_lambda;
  
  A      = 10 * (A_0 + sex *      A_sex +      A_i);
  mu     =      mu_0 + sex *     mu_sex +     mu_i;
  lambda =  lambda_0 + sex * lambda_sex + lambda_i;
}

model {
  vector[M] x;

  for (i in 1:(M)){
      x[i] = logistic(time[i], A[ind[i]], mu[ind[i]], lambda[ind[i]]);
  }
  if(run_estimation==1){
    y ~ normal(x, sigma);
  }
  
  mu_0 ~ normal(0.0, 0.5);
  A_0 ~ normal(0, 1);
  lambda_0 ~ normal(0.0, 2);
  
  A_sex ~ normal(0.0, 3);
  mu_sex ~ normal(0.0, 0.1);
  lambda_sex ~ normal(0.0, 2);
  
  sigma ~ normal(0, 1);
  
  sigma_A ~ normal(0, 0.5);
  sigma_mu ~ normal(0, 0.1);
  sigma_lambda ~ normal(0, 1);

  A_tilde ~ normal(0, 1);
  mu_tilde ~ normal(0, 1);
  lambda_tilde ~ normal(0, 1);
}

generated quantities {
  vector[M] y_sim;
  vector[M] x_sim;

  for (i in 1:(M)){
      x_sim[i] = logistic(time[i], A[ind[i]], mu[ind[i]], lambda[ind[i]]);
  }
  for(i in 1:M) {
    y_sim[i] = normal_rng(x_sim[i], sigma);
  }
}
