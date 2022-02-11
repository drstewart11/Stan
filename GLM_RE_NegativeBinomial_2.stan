data {
  int n;
  int<lower=0> y[n];
  int<lower=0> p;
  int<lower=0> nyrs;
  int<lower=0> ngrid;
  int<lower=1,upper=ngrid> grid_id[n];
  int<lower=1,upper=nyrs> year_id[n];
  int<lower=0, upper=1> PRIOR_ONLY;
  vector[n] salinity;
  vector[n] doxy;
  vector[n] season;
  vector[n] wtemp;
  vector[n] gear;
}

parameters {
  real beta[p];
  vector[nyrs] gamma;
  vector[ngrid] delta;
  real<lower=0> sigma_gamma;
  real<lower=0> sigma_delta;
  real<lower=0> phi;
}

model {
 target += normal_lpdf(beta[1]|0,0.5);
 target += normal_lpdf(beta[2]|0,0.5);
 target += normal_lpdf(beta[3]|0,0.5);
 target += normal_lpdf(beta[4]|0,0.5);
 target += normal_lpdf(beta[5]|0,0.5);
 target += normal_lpdf(beta[6]|0,0.5);
 target += normal_lpdf(gamma|0,sigma_gamma);
 target += cauchy_lpdf(sigma_gamma|0,3.5);
 target += normal_lpdf(delta|0,sigma_delta);
 target += cauchy_lpdf(sigma_delta|0,3.5);
 target += cauchy_lpdf(phi|0,3);

 if(PRIOR_ONLY ==0){
  for(i in 1:n){
   target += neg_binomial_2_log_lpmf(y[i]|beta[1] + gamma[year_id[n]] + delta[grid_id[n]] + beta[2]*salinity[i] + beta[3]*doxy[i] + beta[4]*season[i] + beta[5]*wtemp[i] + beta[6]*gear[i], phi);
}
}
}

generated quantities{
real ypred[n];
real log_lik[n];

for(i in 1:n){
ypred[i] = neg_binomial_2_log_rng(beta[1] + gamma[year_id[n]] + delta[grid_id[n]] + beta[2]*salinity[i] + beta[3]*doxy[i] + beta[4]*season[i] + beta[5]*wtemp[i] + beta[6]*gear[i], phi);
log_lik[i] = neg_binomial_2_log_lpmf(y[i]|beta[1] + gamma[year_id[n]] + delta[grid_id[n]] + beta[2]*salinity[i] + beta[3]*doxy[i] + beta[4]*season[i] + beta[5]*wtemp[i] + beta[6]*gear[i], phi);
}
}


