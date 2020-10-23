data {
  int  N;
  int  J;
  int  M;
  matrix[N, J] X;
  matrix[M, M] phi0;
}

parameters {
  matrix[M, M] phi; // latent variable covaraince matrix
  matrix[N, M] ksi; //latent variable values
  matrix[J, M] lambda; //factor loadings matrix
  real tau[J]; //intercepts
  real<lower=0> psi[J]; //residual variance
}

model {
  // likelihood for data
  for(i in 1:N){
    X[i, 1] ~ normal(tau[1] + ksi[i,1]*lambda[1], psi[1]);
    X[i, 2] ~ normal(tau[2] + ksi[i,1]*lambda[2], psi[2]);
    X[i, 3] ~ normal(tau[3] + ksi[i,1]*lambda[3], psi[3]);
    X[i, 4] ~ normal(tau[4] + ksi[i,2]*lambda[4], psi[4]);
    X[i, 5] ~ normal(tau[5] + ksi[i,2]*lambda[5], psi[5]);
  }
  // latent variable variance matrix
  for(i in 1:N){
    // prior for ksi
    ksi[i] ~ multi_normal(rep_vector(0, M), phi);
  }
  phi  ~ inv_wishart(2, phi0);
  // prior for measurement model parameters
  tau ~ normal(3, 10);
  psi ~ inv_gamma(5, 10);
  lambda[1] ~ normal(1, .001);
  lambda[2] ~ normal(1, 10);
  lambda[3] ~ normal(1, 10);
  lambda[4] ~ normal(1, .001);
  lambda[5] ~ normal(1, 10);
}
