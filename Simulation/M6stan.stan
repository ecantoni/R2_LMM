data {
  int<lower=0> n;
  int<lower=0> J;
  real y[n];
  real x[n];
  real u[J];
  int county[n];
}
parameters { 
  real a[J];
  real b[J];
  
  real gamma0;
  real gamma1;
  real delta0;
  real delta1;
  
  real<lower=0> sigma_y;
  real<lower=0> sigma_a;
  real<lower=0> sigma_b;
  real rho;
}
transformed parameters {
  real y_hat[n];
  vector[2] g[J]; 
  vector[2] g_hat[J];
  matrix[2,2] sigma_g;
  matrix[2,2] chol_sigma_g;
  
  for (i in 1:n){
    y_hat[i] <- a[county[i]] + b[county[i]]*x[i];
	}

  for (j in 1:J){
    g[j,1] <- a[j];
	g[j,2] <- b[j];
    g_hat[j,1] <- gamma0 + gamma1*u[j];
	g_hat[j,2] <- delta0 + delta1*u[j];
	}
	
	sigma_g[1,1] <- sigma_a*sigma_a;
	sigma_g[1,2] <- rho*sigma_a*sigma_b;
	sigma_g[2,1] <- rho*sigma_a*sigma_b;
	sigma_g[2,2] <- sigma_b*sigma_b;

	chol_sigma_g <- cholesky_decompose(sigma_g);
}
model {
  gamma0 ~ normal(0, 100);
  gamma1 ~ normal(0, 100);
  delta0 ~ normal(0, 100);
  delta1 ~ normal(0, 100);
  
  sigma_y ~ uniform(0, 1000);
  sigma_a ~ uniform(0, 100);
  sigma_b ~ uniform(0, 100);
  rho ~ uniform(-1, 1);
  
  for (i in 1:n){
    y[i] ~ normal(y_hat[i], sigma_y);
	}
	
  for (j in 1:J){
    g[j] ~ multi_normal_cholesky(g_hat[j], chol_sigma_g);
	}
}
generated quantities{
  real e_y[n];
  real e_a[J];
  real e_b[J];	
  real dev;
  
  for (i in 1:n){
	e_y[i] <- y[i] - y_hat[i];
	}

  for (j in 1:J){
    e_a[j] <- a[j] - g_hat[j,1];
	e_b[j] <- b[j] - g_hat[j,2];
	}
	
  dev <- 0;
  for (i in 1:n){
    dev <- dev + (-2)*normal_log(y[i], y_hat[i], sigma_y);
	}
}