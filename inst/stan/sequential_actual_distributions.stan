// A version of the sequential model where the actual probability distributions for each calibrated date are passed to the stan programme
functions {
  real calibrated_age_lpdf(real age, data vector ages, data vector probabilities) {
    //takes:
    //age: an age we want to know the probability of
    //ages: a vector of ages in our distribution
    //probabilities: a vector, the same size of ages, which gives the probability of each age
    int lower = 1;
    int upper = num_elements(ages);
    int index = lower + ((upper-lower) / 2);
    int iteration = 0;
    
    real used_age;
    real used_prob;
    real interp_fraction;
    
    if(num_elements(ages) != num_elements(probabilities)) reject("Number of elements in ages and probabilities unequal in calibrated_age_lpdf");
    
    // First use a binary search to find indices for age
    if(ages[lower] > age) {
      upper = lower; //won't go below lower bound of the calibration curve to avoid -Inf log probabilities
    } else if (ages[upper] < age) {
      lower = upper; //won't go above lower bound of the calibration curve to avoid -Inf log probabilities
    } else {
      while(1==1) {
        if(ages[index] == age) {
          lower = index;
          upper = index;
        } else if (ages[index] < age) {
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
      used_age = ages[lower];
      used_prob = probabilities[lower];
    } else {
      interp_fraction = (age - ages[lower]) / (ages[upper] - ages[lower]);
      used_age = (1 - interp_fraction) * ages[lower] + interp_fraction * ages[upper];
      used_prob = (1 - interp_fraction) * probabilities[lower] + interp_fraction * probabilities[upper];
    }
    

    return log(used_prob);
  }
}
data {
  int<lower=1> num_dates;
  int<lower=1> num_probs;
  vector[num_probs] ages[num_dates];
  vector[num_probs] probs[num_dates];
}
transformed data {
  real max_age = 0;
  for(i in 1:num_dates) {
    if (max_age<max(ages[i])) max_age=max(ages[i]);
  }
}
parameters {
  ordered[num_dates] true_ages;
}
model {
  true_ages ~ uniform(0, max_age);
  for(i in 1:num_dates) {
    true_ages[i] ~ calibrated_age(ages[i], probs[i]);
  }
}

