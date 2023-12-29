data {
  int<lower=1> N; // the number of dates
  vector[N] dates;
  vector<lower=0>[N] error;
}
parameters {
  ordered[N] true_ages;
  vector<lower=0>[N] true_error;
}
model {
  true_error ~ normal(error, 1);
  dates ~ normal(true_ages, true_error);
}
