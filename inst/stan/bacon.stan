data {
  int<lower=1> numdets; //j, the number of radiocarbon determinations
  vector[numdets] det_depths; // y, the actual depth measurements
  vector[numdets] det_ages;   // the means of the calibrated ages of the measurements
  vector[numdets] det_sd;     // the sds of the calibrated ages of the measurements
  
  int<lower=1> nummodelledrates; // K, the number of rates to be modelled
  vector[nummodelledrates+1] model_depths;         // the depths to be modelled
  real depth_difference; // delta c
  int<lower=0,upper=1> use_normal;
  
  real a_alpha;  //parameters for gamma distribution
  real b_alpha;
  real a_w;      //parameters for beta distribution
  real b_w;
  
  int<lower=0> t_df; // degrees of freedom for t distribution
}
parameters {
  real<lower=0, upper=1> memory; // w
  vector<lower=0>[nummodelledrates] alpha; //noise on rate AR series
  real<lower=0> min_age;
}
transformed parameters {
  vector<lower=0>[nummodelledrates] rates; // x
  vector<lower=0>[nummodelledrates+1] modelled_ages; // ages at each c. - there's one more than there are rates because there's a c at each end.
  vector<lower=0>[numdets] interpolated_ages; // ages for each applicable c interpolated to values of y and d

  rates[1] = alpha[1];
  for(i in 2:nummodelledrates) {
    rates[i] = memory * rates[i-1] + (1-memory) * alpha[i];
  }
  
  modelled_ages[1] = min_age;
  for(i in 2:nummodelledrates+1) {
    modelled_ages[i] = min_age;
    for(j in 1:(i-1)) {
      modelled_ages[i] += rates[j] * depth_difference;
    }
  }

  for(i in 1:numdets) {
    interpolated_ages[i] = min_age;
    for(j in 2:nummodelledrates +1 ) {
      if(model_depths[j] <= det_depths[i]) {
        interpolated_ages[i] += depth_difference * rates[j-1];
      } else {
        interpolated_ages[i] += (det_depths[i] - model_depths[j-1]) * rates[j-1];
        break;
      }
    }
  }
  

}

model {
  if(use_normal==1) {
    det_ages[1:numdets] ~ normal(interpolated_ages[1:numdets], det_sd[1:numdets]);
  } else {
    det_ages[1:numdets] ~ student_t(t_df, interpolated_ages[1:numdets], det_sd[1:numdets]);
  }
  min_age ~ normal(det_ages[1], det_sd[1]);
  memory ~ beta(a_w, b_w);
  alpha[1:nummodelledrates] ~ gamma(a_alpha, b_alpha);
}

