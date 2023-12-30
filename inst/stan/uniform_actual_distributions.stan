// A version of the uniform model where the actual probability distributions for each calibrated date are passed to the stan programme
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

  real calibrated_age2_lpdf(real age, matrix calib_ages) {
    return log(binary_search(age, calib_ages[,1], calib_ages[,2]));
  }
  
  real calibrated_age_lpdf(real age, data vector ages, data vector probabilities) {
    //takes:
    //age: an age we want to know the probability of
    //ages: a vector of ages in our distribution
    //probabilities: a vector, the same size of ages, which gives the probability of each age
    return log(binary_search(age, ages, probabilities));
  }
}
data {
  int<lower=1> num_dates;
  int<lower=1> num_probs;
  vector[num_probs] ages[num_dates];
  vector[num_probs] probs[num_dates];
  vector[num_dates] depth;
  //the calibration curve
  int<lower=1> calcurve_size;
  vector[calcurve_size] calibrated_ages;
  vector[calcurve_size] c14_ages;
  vector[calcurve_size] cc_sd;
}
parameters {
  real rate;
//  vector[num_dates] det_ages_dist;
  real first_age;
}
transformed parameters {
  ordered[num_dates] modelled_age;

//  modelled_age[1] = det_ages_dist[1];
  modelled_age[1] = first_age;
  for(i in 2:num_dates) {
    modelled_age[i] = modelled_age[i-1] + rate * (depth[i] - depth [i-1]);
  }
}
model {
  first_age ~ calibrated_age(ages[1], probs[1]);
  for(i in 1:num_dates) {
    modelled_age[i] ~ calibrated_age(ages[i], probs[i]);
  }
}

