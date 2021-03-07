data {
  int<lower=0> n;
  real y[n];
  real x[n];
}
parameters {
  real a;
  real b;
  real<lower=0> sigma_y;
}
transformed parameters {
  real y_hat[n];
  
  for (i in 1:n){
    y_hat[i] <- a + b*x[i];
	}
}
model {
  a ~ normal(0, 100);
  b ~ normal(0, 100);
  sigma_y ~ uniform(0, 1000);
  
  for (i in 1:n){
    y[i] ~ normal(y_hat[i], sigma_y);
	}
}
generated quantities{
  real e_y[n];
  real dev;
  
  for (i in 1:n){
	e_y[i] <- y[i] - y_hat[i];
	}

  dev <- 0;
  for (i in 1:n){
    dev <- dev + (-2)*normal_log(y[i], y_hat[i], sigma_y);
	}
}