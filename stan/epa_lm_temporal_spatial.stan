data {
  int<lower=0> n;
  int<lower=0> p;
  int<lower=0> t;
  int<lower=0> c;  
  matrix[n, p] X;
  matrix[n, t] Z1;
  matrix[n, c] Z2;  
  matrix<lower = 0, upper = 1>[c, c] W;
  real y[n];
}
transformed data{
  vector[c] zeros;
  matrix<lower = 0>[c,c] D;
  {
    vector[c] W_rowsums;
    for (i in 1:c) {
      W_rowsums[i] = sum(W[i, ]);
    }
    D = diag_matrix(W_rowsums);
  }
  zeros = rep_vector(0, c);
}
parameters {
  vector[p] beta;
  vector[t] omega;
  vector[c] theta;
  real<lower=-0.999,upper=0.999> phi;
  real<lower = 0, upper = 1> rho;
  real<lower=0> sigma;
  real<lower=0> sigma_w;
  real<lower=0> tau;    
}
transformed parameters {
  vector[n] mu;
  vector[t] eta;
  eta[1] = 0;
  for(i in 2:t) {
	eta[i] = omega[i] + phi*eta[i-1];
  }  
  mu = X * beta + Z1 * eta + Z2 * theta;
}

model {
  //Hierarchical Model  
  omega ~ normal(0, sigma_w);
  theta ~ multi_normal_prec(zeros, tau * (D - rho * W)); 
  
  //Priors
  beta ~ normal(0,2);  
  sigma ~ cauchy(0, 5);
  phi ~ normal(0,1);
  sigma_w ~ cauchy(0, 5);  
  tau ~ gamma(2, 2);
  
  //Data Model
  log(y) ~ normal(mu, sigma);
}

