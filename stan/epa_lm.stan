data {
  int<lower=0> n;
  int<lower=0> p;
  matrix[n, p] X;
  vector[n] y;
}
parameters {
  vector[p] beta;  
  real<lower=0> sigma;
}
transformed parameters {  
  vector[n] mu;
  mu = X * beta;
}
model {
  //Hierachical Model
  beta ~ normal(0,2);
  
  //Prior
  sigma ~ cauchy(0, 1);
  
  //Data Model
  log(y) ~ normal(mu, sigma);
}
