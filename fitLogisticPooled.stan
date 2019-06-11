functions {
    // logistic function
    real logistic(real time, real A, real mu, real lambda) {
        real x;
        x = A/(1 + exp(4 * mu * (lambda - time)/A + 2));
        return x;
    }
}

data {
  int<lower=0> M; // numero de observacoes
  vector[M] y;    // pesos
  vector[M] sex;  // sexo das observacoes
  vector[M] time; // dia das observacoes
}

parameters {
  real<lower=0> mu;
  real<lower=0>  A;
  real lambda;
  
  real mu_s;
  real A_s;
  real lambda_s;

  real<lower=0> sigma;
}

model {
  vector[M] x;
  for(m in 1:M){
    x[m] = logistic(time[m], 10 * (A + A_s * sex[m]), mu + mu_s * sex[m], lambda + lambda_s * sex[m]);
  }
  y ~ normal(x, sigma);
  
  mu ~ normal(0, 1);
  A ~ normal(0, 2);
  lambda ~ normal(0, 2);

  mu_s ~ normal(0, 0.5);
  A_s ~ normal(0, 0.5);
  lambda_s ~ normal(0, 0.5);
  
  sigma ~ normal(0, 1);
}





