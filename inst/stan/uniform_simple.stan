data {
  int N; //number of measurements
  vector[N] depth;
  vector[N] age;
  vector[N] age_sd;
}
parameters {
  real rate;
  real first_age;
}
transformed parameters {
  ordered[N] modelled_age;
  modelled_age[1] = first_age;
  for(i in 2:N) {
    modelled_age[i] = modelled_age[i-1] + rate * (depth[i] - depth [i-1]);
  }
}
model {
  first_age ~ normal(age[1], age_sd[1]);
  age ~ normal(modelled_age, age_sd);
}

// data {
//   int N; //number of measurements 
//   vector[N] depth;
//   vector[N] age;
//   vector[N] age_sd;
// }
// parameters {
//   ordered[N] post_age;
// //  vector[N] post_age;
//   real slope;
//   real intercept;
//   real<lower=0> sigma;
// }
// model {
//   post_age ~ normal(age, age_sd);
//   post_age ~ normal(intercept + slope * depth, sigma);
// }



// //
// // A simple uniform model
// // Takes two radiocarbon dates (mean & sd)
// // and two depths
// // And calculates a straight line
// 
// data {
//   real depth_top;
//   real depth_bottom;
//   real age_top;
//   real age_bottom;
//   real age_sd_top;
//   real age_sd_bottom;
// }
// transformed data {
//   real approx_intercept;
//   real approx_slope;
//   
//   approx_slope = (age_bottom - age_top) / (depth_bottom - depth_top);
//   approx_intercept = age_bottom - depth_bottom * (approx_slope);
// }
// parameters {
//   ordered[2] post_ages;
//   real slope;
//   real intercept;
//   real<lower=0> sigma;
// }
// model {
//   intercept    ~ normal(approx_intercept, 10 * approx_intercept);
//   slope        ~ normal(approx_slope, 10 * approx_slope);
//   age_top ~ normal(post_ages[1], age_sd_top);
//   age_bottom ~ normal(post_ages[2], age_sd_bottom);
//   post_ages[1] ~ normal(intercept + slope * depth_top, sigma);
//   post_ages[2] ~ normal(intercept + slope * depth_bottom, sigma);
// }
// 
// 


