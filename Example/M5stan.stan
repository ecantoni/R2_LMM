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
}
transformed parameters {
  real y_hat[n];
  real a_hat[J];
  real b_hat[J];
  
  for (i in 1:n){
    y_hat[i] <- a[county[i]] + b[county[i]]*x[i];
	}

  for (j in 1:J){
    a_hat[j] <- gamma0 + gamma1*u[j];
	b_hat[j] <- delta0 + delta1*u[j];
	}
}
model {
  gamma0 ~ normal(0, 100);
  gamma1 ~ normal(0, 100);
  delta0 ~ normal(0, 100);
  delta1 ~ normal(0, 100);
  
  sigma_y ~ uniform(0, 1000);
  sigma_a ~ uniform(0, 100);
  sigma_b ~ uniform(0, 100);
  
  for (i in 1:n){
    y[i] ~ normal(y_hat[i], sigma_y);
	}
	
  for (j in 1:J){
    a[j] ~ normal(a_hat[j], sigma_a);
    b[j] ~ normal(b_hat[j], sigma_b);
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
    e_a[j] <- a[j] - a_hat[j];
	e_b[j] <- b[j] - b_hat[j];
	}

  dev <- 0;
  for (i in 1:n){
    dev <- dev + (-2)*normal_log(y[i], y_hat[i], sigma_y);
	}
}