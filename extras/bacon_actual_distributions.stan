functions {
  real binary_search(real search, vector index_values, vector output) {
    //searches for values in output corresponding to the values in the index identified by the search parameter
    //if no exact match is found, interpolates the two output values
    int lower = 1;
    int upper = num_elements(index_values);
    int index = lower + ((upper - lower) / 2);
    int iteration = 0;
    
    real used_index_value;
    real used_output;
    real interp_fraction;
    
    if(num_elements(index_values) != num_elements(output)) reject("Number of elements in index and output vectors unequal in binary_search");
    
    if(index_values[lower] > search) {
      upper = lower; //won't go below lower bound of the calibration curve to avoid -Inf log probabilities
    } else if (index_values[upper] < search) {
      lower = upper; //won't go above lower bound of the calibration curve to avoid -Inf log probabilities
    } else {
      while(1==1) {
        if(index_values[index] == search) {
          lower = index;
          upper = index;
        } else if (index_values[index] < search) {
          lower = index;
        } else {
          upper = index;
        }
        if (upper - lower <= 1 || iteration > 50) break;
        index = lower + ((upper-lower) / 2);
        iteration = iteration + 1;
      }
    }

    if(upper == lower) {
      used_index_value = index_values[lower];
      used_output = output[lower];
    } else {
      interp_fraction = (search - index_values[lower]) / (index_values[upper] - index_values[lower]);
      used_index_value = (1 - interp_fraction) * index_values[lower] + interp_fraction * index_values[upper];
      used_output = (1 - interp_fraction) * output[lower] + interp_fraction * output[upper];
    }
    
    return used_output;
    
  }

  matrix calibrate(real measurement, real measurement_sd, vector calibrated_ages, vector c14_ages, vector cc_sd) {
    // calibrate a single C14 measurement
    // returns an matrix containing ages in column 1 and probabilities in column 2
    matrix[num_elements(calibrated_ages), 2] output;
    real total_prob;

    if(num_elements(calibrated_ages) != num_elements(c14_ages) ||
      num_elements(calibrated_ages) != num_elements(cc_sd) ||
      num_elements(c14_ages) != num_elements(cc_sd))
        reject("Number of elements in index and output vectors unequal in calibrate");
    
    total_prob = 0;
    
    for(i in 1:num_elements(calibrated_ages)) {
      output[i, 1] = calibrated_ages[i];
      output[i, 2] = exp(normal_lpdf(c14_ages[i] | measurement, sqrt(measurement_sd^2 + cc_sd[i]^2)));
//      total_prob = total_prob + output[i, 2];
    } 
  
    // total_prob = sum(output[, 2]);
    // 
    // for(i in 1:num_elements(output[,2])) {
    //   output[i, 2] = output[i, 2] / total_prob;
    // }

    return output;
  }
  
  real calibrated_age_lpdf(real age, vector ages, vector probabilities) {
    return log(binary_search(age, ages, probabilities));
  }
  
  real calibrated_age2_lpdf(real age, matrix calib_ages) {
    return log(binary_search(age, calib_ages[,1], calib_ages[,2]));
  }
}
data {
  int<lower=1> numdets; //j, the number of radiocarbon determinations
  vector[numdets] det_depths; // y, the actual depth measurements
  //input values are the full distributions of the calibrated input dates
  int<lower=1> num_probs; // the number of probability values stored for each depth
  vector[num_probs] ages[numdets];
  vector[num_probs] probs[numdets];
  vector[numdets] det_meas_age;
  vector[numdets] det_meas_sd;
  vector[numdets] det_calib_age;
  vector[numdets] det_calib_sd;
  vector[numdets] min_ages;
  vector[numdets] max_ages;
  vector[numdets] min_c14_ages;
  vector[numdets] max_c14_ages;
  real abs_min_age;
  real abs_max_age;
  real abs_min_c14_age;
  real abs_max_c14_age;

  int<lower=1> nummodelledrates; // K, the number of rates to be modelled
  vector[nummodelledrates+1] model_depths;         // the depths to be modelled
  real depth_difference; // delta c
  int<lower=0,upper=1> use_normal;
  
  real a_alpha;  //parameters for gamma distribution
  real b_alpha;
  real a_w;      //parameters for beta distribution
  real b_w;
  
  int<lower=0> t_df; // degrees of freedom for t distribution
  
  //the calibration curve
  int<lower=1> calcurve_size;
  vector[calcurve_size] calibrated_ages;
  vector[calcurve_size] c14_ages;
  vector[calcurve_size] cc_sd;

}

parameters {
  real<lower=0, upper=1> memory; // w
  vector<lower=0>[nummodelledrates] alpha; //noise on rate AR series
  real<lower=min_c14_ages[1], upper=max_c14_ages[1]> min_age;
  vector<lower=abs_min_age, upper=abs_max_age>[numdets] det_ages_dist;
  vector<lower=abs_min_age, upper=abs_max_age>[numdets] interpolated_ages_calib;
}
transformed parameters {
  vector<lower=0>[nummodelledrates] rates; // x
  vector[nummodelledrates+1] modelled_ages; // ages at each c. - there's one more than there are rates because there's a c at each end.
  vector[numdets] interpolated_ages; // ages for each applicable c interpolated to values of y and d

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
  for(i in 1:numdets) {
    interpolated_ages_calib[i] ~ calibrated_age2(calibrate(interpolated_ages[i], det_meas_sd[i], calibrated_ages, c14_ages, cc_sd));
    det_ages_dist[i] ~ calibrated_age(ages[i], probs[i]);
    if(use_normal==1) {
      det_ages_dist[i] ~ normal(interpolated_ages_calib[i], det_calib_sd[i]);
    } else {
      det_ages_dist[i] ~ student_t(t_df, interpolated_ages_calib[i], 2 * det_calib_sd[i]);
    }
  }
  min_age ~ normal(det_meas_age[1], det_meas_sd[1]); 
  memory ~ beta(a_w, b_w);
  alpha[1:nummodelledrates] ~ gamma(a_alpha, b_alpha);
}


