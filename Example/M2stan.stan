data {
  int<lower=0> n;
  int<lower=0> J;
  real y[n];
  int county[n];
}
parameters {
  real a[J];
  real mu_a;
  
  real<lower=0> sigma_y;
  real<lower=0> sigma_a;
}
transformed parameters {
  real y_hat[n];
  
  for (i in 1:n){
    y_hat[i] <- a[county[i]];
	}
}
model {
  a ~ normal(0, 100);
  mu_a ~ normal(0, 100);
  
  sigma_y ~ uniform(0, 1000);
  sigma_a ~ uniform(0, 100);
  
  for (i in 1:n){
    y[i] ~ normal(y_hat[i], sigma_y);
	}
	
  for (j in 1:J){
    a[j] ~ normal(mu_a, sigma_a);
	}
}
generated quantities{
  real e_y[n];	
  real e_a[J];
  real dev;
  
  for (i in 1:n){
	e_y[i] <- y[i] - y_hat[i];
	}

  for (j in 1:J){
    e_a[j] <- a[j] - mu_a;
	}

  dev <- 0;
  for (i in 1:n){
    dev <- dev + (-2)*normal_log(y[i], y_hat[i], sigma_y);
	}
}