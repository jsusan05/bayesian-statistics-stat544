data {
  int<lower=0> n;
  int<lower=0> p;
  int<lower=0> t;
  matrix[n, p] X;
  matrix[n, t] Z1;
  real y[n];  
}
parameters {
  vector[p] beta;
  vector[t] omega;
  real<lower=-0.999,upper=0.999> phi;
  real<lower=0> sigma;
  real<lower=0> sigma_w;    
}
transformed parameters {
  vector[n] mu;
  vector[t] eta;
  eta[1] = 0;
  for(i in 2:t) {
    eta[i] = omega[i] + phi*eta[i-1];
  }
  mu = X * beta + Z1 * eta;
}
model {   
  //Hierarchical Model
  beta ~ normal(0,2);
  omega ~ normal(0, sigma_w);
  phi ~ normal(0,1); 
  
  //Priors
  sigma ~ cauchy(0, 5);
  sigma_w ~ cauchy(0, 5); 

  //Data Model
  log(y) ~ normal(mu, sigma);
}
