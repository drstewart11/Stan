data {
  int n;
  int<lower=0> y[n];
  int<lower=0> p;
  int nyrs;
  int g[n];
  int<lower=0, upper=1> PRIOR_ONLY;
  vector[n] salinity;
  vector[n] doxy;
  vector[n] season;
  vector[n] wtemp;
  vector[n] wdepth;
}

parameters {
  real beta[p];
  real a[nyrs];
  real<lower=0,upper=10> sigma;
}

model {
 target += normal_lpdf(beta[1]|0,0.5);
 target += normal_lpdf(beta[2]|0,0.5);
 target += normal_lpdf(beta[3]|0,0.5);
 target += normal_lpdf(beta[4]|0,0.5);
 target += normal_lpdf(beta[5]|0,0.5);
 target += normal_lpdf(beta[6]|0,0.5);
 target += normal_lpdf(a|0,sigma);

 if(PRIOR_ONLY ==0){
  for(i in 1:n){
   target += poisson_lpmf(y[i]|beta[1] + a[g[n]] + beta[2]*salinity[i] + beta[3]*doxy[i] + beta[4]*season[i] + beta[5]*wtemp[i] + beta[6]*wdepth[i]);
}
}
}

generated quantities{
real ypred[n];
real log_lik[n];

for(i in 1:n){
ypred[i] = poisson_rng(beta[1] + a[g[n]] + beta[2]*salinity[i] + beta[3]*doxy[i] + beta[4]*season[i] + beta[5]*wtemp[i] + beta[6]*wdepth[i]);
log_lik[i] = poisson_lpmf(y[i]|beta[1] + a[g[n]] + beta[2]*salinity[i] + beta[3]*doxy[i] + beta[4]*season[i] + beta[5]*wtemp[i] + beta[6]*wdepth[i]);
}
}


