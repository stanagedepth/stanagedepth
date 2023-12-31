library(rstan)
library(dplyr)
library(tidyr)
library(stringr)
library(posterior)
library(rintcal)


#' @import rstan
#' @import dplyr
#' @import tidyr
#' @import stringr
#' @import posterior
#' @import rintcal


#' @title Some sample age-depth data
#'
#' @description `sample_data` returns a data frame of synthetic age-depth data which can be used as an input to the run_model function
#' @export
sample_data <- function() {
  data.frame(
    label  = c('top', 'd2', 'd3', 'd4', 'd5','d6','d7', 'd8', 'd9', 'bottom'),
    age = c(395, 485, 613, 595, 743, 809, 848, 887, 990, 1310),
    error = c(25, 30, 25, 25, 25, 25, 25, 25, 25, 30),
    depth = c(0.20, 0.26, 0.32, 0.38, 0.44, 0.49, 0.53, 0.57, 0.7, 0.8)
  )
}

calibration_curves <- function() {
  list(
    IntCal20 = 1
  )
}

get_calibrated_date <- function(uncalibrated_age, error, curve = 'IntCal20') {
  calibrated_distribution <- as.data.frame(rintcal::caldist(
    uncalibrated_age,
    error,
    calibration_curves()[[curve]],
    cc.resample = 1
  ))

  names(calibrated_distribution) <- c('cal_bp', 'prob')

  mean_calibrated_age <- sum(calibrated_distribution$cal_bp * calibrated_distribution$prob)

  squared_deviations <- (calibrated_distribution$cal_bp - mean_calibrated_age)^2
  weighted_squared_deviations <- squared_deviations * calibrated_distribution$prob
  sum_ws_deviations <- sum(weighted_squared_deviations)
  sd_calibrated_age <- sqrt(sum_ws_deviations)

  min_calibrated_age <- min(calibrated_distribution[, 1])
  max_calibrated_age <- max(calibrated_distribution[, 1])


  calibrated_date <- list(
    calibrated_distribution = calibrated_distribution,
    uncalibrated_age = uncalibrated_age,
    age_error = error,
    mean_calibrated_age = mean_calibrated_age,
    sd_calibrated_age = sd_calibrated_age,
    min_calibrated_age = min_calibrated_age,
    max_calibrated_age = max_calibrated_age
  )


  class(calibrated_date) <- 'calibdate'

  calibrated_date
}


get_date_samples <- function(x, size = 1000) {
  if (class(x) != 'calibdate') stop('x in get_date_samples not of proper type')
  sample(x$calibrated_distribution$cal_bp,
         size = size,
         replace = TRUE,
         prob = x$calibrated_distribution$prob)
}


preprocess_input <- function(x, curve = 'IntCal20', samples = 1000) {
  # Takes an input df containing: label, c14 age, c14 age error, depth;
  # and adds calibrated dates
  # Names of data frame to be label, age, error, depth
  # Added columns are calib_date_mean, calib_sd, calibrated_dist
  if (any(names(x) != c('label', 'age', 'error', 'depth'))) stop('Names of input not as expected')
  calib_date_mean <- numeric(nrow(x))
  calib_date_sd   <- numeric(nrow(x))
  calibrated_dist <- vector('list', nrow(x))
  calibrated_dist_samples <- vector('list', nrow(x))
  for (i in seq_len(nrow(x))) {
    cdate <- get_calibrated_date(x$age[i], x$error[i], curve)
    calib_date_mean[i] <- cdate$mean_calibrated_age
    calib_date_sd[i] <- cdate$sd_calibrated_age
    calibrated_dist[i] <- list(cdate$calibrated_distribution)
    calibrated_dist_samples[i] <- list(get_date_samples(cdate, samples))
  }
  cbind(x, calib_date_mean, calib_date_sd, I(calibrated_dist), I(calibrated_dist_samples))
}


get_bacon_rates <- function(stan_output, model_depths) {

  rates <- posterior::as_draws_df(stan_output) %>%
    mutate(iteration = row_number()) %>%
    select(starts_with('rates')|'iteration') %>%
    pivot_longer(starts_with('rates')) %>%
    rename(depth_index = name, rate = value) %>%
    mutate(depth_index = as.numeric(str_replace_all(str_replace_all(depth_index, 'rates\\[', ''), ']', ''))) %>%
    inner_join(model_depths, by = 'depth_index') %>%
    arrange(iteration, depth_index) %>%
    group_by(iteration) %>%
    mutate(cum_rate = cumsum(rate)) %>%
    ungroup() %>%
    mutate(cum_rate = cum_rate - rate)

  return(rates)
}

get_bacon_ages <- function(x, stan_output, rates, slice_thick, output_points, model_depths) {

  min_ages <- pull(posterior::as_draws_df(stan_output), 'min_age')

  ages_at_depth <- function(d) {
    if(d == max(model_depths$depth)) {
      max_di = max(rates$depth_index)
      top_rate <- rates %>% filter(depth_index == max_di)
       min_ages + (slice_thick * (pull(top_rate, cum_rate) + pull(top_rate, rate)))
    } else {
      top_rate <- rates %>% filter(depth_index == max(which(model_depths$depth <= d)))
      min_ages + (slice_thick * pull(top_rate, cum_rate)) + ((d - pull(top_rate, depth)) * pull(top_rate, rate))
    }
  }

  output_depths <- seq(from = min(x$depth), to = max(x$depth), length.out = output_points)
  output_means <- numeric(output_points)
  output_mins <- numeric(output_points)
  output_maxs <- numeric(output_points)
  output_sds  <- numeric(output_points)


  for(i in seq_along(output_depths)) {
    ages <- ages_at_depth(output_depths[i])
    output_means[i] <- mean(ages)
    hdi <- HDInterval::hdi(ages, credMass = 0.95)
    output_mins[i] <- hdi[1]
    output_maxs[i] <- hdi[2]
    output_sds[i] <- sd(ages)
  }

  as.data.frame(list(depth = output_depths, min = output_mins, mean = output_means, max = output_maxs, sd = output_sds))
}

bacon_full_reprocess_ages <- function(x, output_ages, output_points) {

  output_depths <- seq(from = min(x$depth), to = max(x$depth), length.out = output_points)
  output_means <- numeric(output_points)
  output_mins <- numeric(output_points)
  output_maxs <- numeric(output_points)
  output_sds  <- numeric(output_points)


  for(i in 1:nrow(output_ages)) {
    output_depths <- output_ages$depth
    print(paste("i",i,"mean",output_ages$mean[i],"sd",output_ages$sd[i]))
    age <- get_calibrated_date(output_ages$mean[i], output_ages$sd[i])
    print(age$mean_calibrated_age)
    output_means[i] <- age$mean_calibrated_age
    hdi <- HDInterval::hdi(
      sample(age$calibrated_distribution$cal_bp, size = 10000, prob = age$calibrated_distribution$prob, replace = TRUE),
      credMass = 0.95)
    output_mins[i] <- hdi[1]
    output_maxs[i] <- hdi[2]
    output_sds[i] <- age$sd_calibrated_age
  }

  as.data.frame(list(depth = output_depths, min = output_mins, mean = output_means, max = output_maxs, sd = output_sds))
}

run_bacon_simple <- function(x, chains = 4, warmup = 1000, iter = 3000, modelled_points = 20, output_points = 100, a_alpha = 1, b_alpha = 2, a_w = 2, b_w = 2, t_df = 3, use_normal = 0) {

  model_depths <- as.data.frame(list(depth_index = seq_len(modelled_points), depth = seq(from = min(x$depth), to = max(x$depth), length.out = modelled_points)))
  slice_thick <- pull(model_depths, depth)[2] - pull(model_depths, depth)[1]

  stan_model_data <- list(
    numdets = nrow(x),
    det_depths = x$depth,
    det_ages = x$calib_date_mean,
    det_sd = x$calib_date_sd,
    nummodelledrates = modelled_points - 1,
    model_depths = model_depths$depth,
    depth_difference = slice_thick,
    a_alpha = a_alpha,
    b_alpha = b_alpha,
    a_w = a_w,
    b_w = b_w,
    t_df = t_df,
    use_normal = use_normal
  )

  model <- stanmodels$bacon

  stan_output <- rstan::sampling(
    object = model,
    data = stan_model_data,
    chains = chains,
    warmup = warmup,
    iter = iter,
    cores = parallel::detectCores())

  rates <- get_bacon_rates(stan_output, model_depths)
  bacon_ages <- get_bacon_ages(x, stan_output, rates, slice_thick, output_points, model_depths)

  age_distrib <- posterior::as_draws_df(stan_output) %>% select(starts_with('interpolated_ages'))
  names(age_distrib) <- x$label
  age_distrib <- age_distrib %>% pivot_longer(everything(), names_to = 'label', values_to = 'age') %>%
    inner_join(select(x, label, depth), by = 'label')

  output <- list(
    input = x,
    modelled_depths = model_depths$depth,
    stan_output = stan_output,
    output_ages = bacon_ages,
    rates = rates,
    age_distrib = age_distrib,
    model_family = 'bacon',
    implementation_family = 'simple',
    implementation = 'bacon_simple'
    )

  class(output) <- append(class(output), 'agedepthmodel')

  return(output)

}


run_bacon <- function(x, chains = 4, warmup = 1000, iter = 3000, modelled_points = 20, output_points = 100, a_alpha=1.5, b_alpha=.075, a_w=5, b_w=5, t_df=3, use_normal = 0) {

  bacon_simple <- run_bacon_simple(x,
                                   chains = 4,
                                   warmup = 500,
                                   iter = 1200,
                                   modelled_points = modelled_points,
                                   a_alpha = a_alpha,
                                   b_alpha = b_alpha,
                                   a_w = a_w,
                                   b_w = b_w,
                                   t_df = t_df,
                                   use_normal = use_normal)

  bacon_simple_df <- tibble(posterior::as_draws_df(bacon_simple$stan_output))

  model_depths <- as.data.frame(list(depth_index = seq_len(modelled_points), depth = seq(from = min(x$depth), to = max(x$depth), length.out = modelled_points)))
  slice_thick <- pull(model_depths, depth)[2] - pull(model_depths, depth)[1]

  distribs <- x$calibrated_dist
  num_dates <- nrow(x)
  num_probs <- max(vapply(distribs, nrow, numeric(1)))

  ages_matrix <- matrix(nrow = num_probs, ncol = num_dates)
  probs_matrix <- matrix(nrow = num_probs, ncol = num_dates)

  for (i in 1:num_dates) {
    dates <- distribs[[i]]$cal_bp
    probs <- distribs[[i]]$prob
    shortfall <- num_probs - length(probs)
    ages_matrix[, i] <- c(rep(-1, shortfall), dates)
    probs_matrix[, i] <- c(rep(0, shortfall), probs)
  }

  cc <- as.data.frame(ccurve(resample=10))
  names(cc) <- c('V1', 'V2', 'V3')

  min_ages = numeric(num_dates)
  max_ages = numeric(num_dates)
  min_c14_ages = numeric(num_dates)
  max_c14_ages = numeric(num_dates)

  scale = c(min = 0.9, max = 1.1)

  for(i in 1:num_dates) {
    interpolated_ages_simple <- bacon_simple_df %>%
      select(starts_with('interpolated_ages')) %>%
      pivot_longer(everything()) %>%
      group_by(name) %>%
      summarise(min_age = min(value),
                max_age = max(value)) %>%
      ungroup() %>%
      mutate(name = as.numeric(str_replace_all(str_replace_all(name, 'interpolated_ages\\[', ''), ']' , ''))) %>%
      arrange(name)
    min_ages[i] <- scale['min'] * min(c(interpolated_ages_simple[i, ]$min_age, min((x[i,])$calibrated_dist_samples[[1]])))
    max_ages[i] <- scale['max'] * max(c(interpolated_ages_simple[i, ]$max_age, max((x[i,])$calibrated_dist_samples[[1]])))
    min_c14_ages[i] <- rintcal::calBP.14C(min_ages[i])[1] - 3* calBP.14C(min_ages[i])[2]
    max_c14_ages[i] <- rintcal::calBP.14C(max_ages[i])[1] + 3* calBP.14C(max_ages[i])[2]
  }

  abs_min_age <- min(min_ages)
  abs_max_age <- max(max_ages)
  abs_min_c14_age <- min(min_c14_ages)
  abs_max_c14_age <- max(max_c14_ages)

   init_values <- function() {

     # Only use the more likely bacon outputs
     bacon_simple_df <- bacon_simple_df %>% slice_min(lp__, prop = 0.5)

     rc14age <- function(x) {
       c14_age <- rintcal::calBP.14C(x)
       rnorm(1, c14_age[1], c14_age[2])
     }
     get_bacon_values <- function(start, row) {
       bacon_simple_df %>% select(starts_with(start)) %>% slice(row) %>% pivot_longer(everything()) %>% pull(value)
     }

     det_ages_iv <- numeric(num_dates)
     for(i in 1:num_dates) {
       det_ages_iv[i] <- rnorm(1, x$calib_date_mean[i], x$calib_date_sd[i])
     }

     row <- floor(runif(1, 1, nrow(bacon_simple_df)))

     interpolated_ages_calib <- get_bacon_values('interpolated_ages', row)
     interpolated_ages <- numeric(length(interpolated_ages_calib))
     interpolated_ages[1] = rc14age(interpolated_ages_calib[1])
     for(i in 2:length(interpolated_ages_calib)) {
       age <- rc14age(interpolated_ages_calib[i])
       count <- 0
       while(age < interpolated_ages[i - 1]) {
         age <- rc14age(interpolated_ages_calib[i])
         count <- count + 1
         if(count>50) {
           age <- interpolated_ages[i - 1]
           break
         }
       }
       interpolated_ages[i] <- age
       interpolated_ages_calib[i] <- get_calibrated_date(age, 0)$mean_calibrated_age
     }

     modelled_ages_calib <- get_bacon_values('modelled_ages', row)
     modelled_ages <- numeric(length(modelled_ages_calib))
     modelled_ages[1] <- rc14age(modelled_ages_calib[1])
     for(i in 2:length(modelled_ages_calib)) {
       age <- rc14age(modelled_ages_calib[i])
       count <- 0
       while(age < modelled_ages[i - 1]) {
         age <- rc14age(modelled_ages_calib[i])
         count <- count + 1
         if(count>50) {
           age <- modelled_ages[i - 1]
           break
         }
       }
       modelled_ages[i] <- age
       modelled_ages_calib[i] <- get_calibrated_date(age, 0)$mean_calibrated_age
     }

     memory <- get_bacon_values('memory', row)
     rates <- get_bacon_values('rates', row)
     alpha <- get_bacon_values('alpha', row)
     min_age <- rc14age(get_bacon_values('min_age', row))


     list(
       min_age = min_age,
       det_ages_dist = det_ages_iv,
       interpolated_ages_calib = interpolated_ages_calib,
       interpolated_ages = interpolated_ages,
       modelled_ages = modelled_ages
     )
   }


  stan_model_data <- list(
    numdets = num_dates,
    det_depths = x$depth,
    num_probs = num_probs,
    ages = t(ages_matrix),
    probs = t(probs_matrix),
    det_calib_age = x$calib_date_mean,
    det_calib_sd = x$calib_date_sd,
    min_ages = min_ages,
    max_ages = max_ages,
    min_c14_ages = min_c14_ages,
    max_c14_ages = max_c14_ages,
    abs_min_age = abs_min_age,
    abs_max_age = abs_max_age,
    abs_min_c14_age = abs_min_c14_age,
    abs_max_c14_age = abs_max_c14_age,
    det_meas_age = x$age,
    det_meas_sd = x$error,
    calcurve_size = nrow(cc),
    calibrated_ages = cc$V1,
    c14_ages = cc$V2,
    cc_sd = cc$V3,
    nummodelledrates = modelled_points - 1,
    model_depths = model_depths$depth,
    depth_difference = slice_thick,
    a_alpha = a_alpha,
    b_alpha = b_alpha,
    a_w = a_w,
    b_w = b_w,
    t_df = t_df,
    use_normal = use_normal
  )


    model_file <- 'extras/bacon_actual_distributions_newstan.stan'
    cmdstanmodel <- cmdstan_model(model_file)
    cmdstan_output <- cmdstanmodel$sample(
      data = stan_model_data,
      chains = chains,
      iter_warmup = warmup,
      iter_sampling = iter - warmup,
      init = init_values,
      parallel_chains = parallel::detectCores(),
      output_dir = 'cmdstanoutput/'
    )
    stan_output <- rstan::read_stan_csv(cmdstan_output$output_files())


  rates <- get_bacon_rates(stan_output, model_depths)
  bacon_ages <- get_bacon_ages(x, stan_output, rates, slice_thick, output_points, model_depths)
  bacon_ages <- bacon_full_reprocess_ages(x, bacon_ages, output_points)


  output <- list(
    input = x,
    modelled_depths = model_depths$depth,
    stan_output = stan_output,
    output_ages = bacon_ages,
    rates = rates,
    age_distrib = 1,
    model_family = 'bacon',
    implementation_family = 'full',
    implementation = 'bacon_full'
  )

  class(output) <- append(class(output), 'agedepthmodel')

  return(output)
}


uniform_postprocess <- function(x, stan_output, implementation_family) {
  stan_data <- posterior::as_draws_df(stan_output)

  output_depths <- x$depth

  N <- 500

  output_ages <- as.data.frame(list(iter = numeric(0), depth = numeric(0), age = numeric(0)))

  for(i in 1:N) {
    for(j in 1:length(output_depths)) {
      row = floor(runif(1, 1, nrow(stan_data)))
      depth = output_depths[j]
      slope = stan_data$rate[row]
      first_age = stan_data$first_age[row]
      age = first_age + slope * (depth - output_depths[1])
      output_ages <- rbind(output_ages, data.frame(iter = i, depth = depth, age = age))
    }
  }

  output_ages_summary <- output_ages %>%
    group_by(depth) %>%
    summarise(mean = mean(age), min = HDInterval::hdi(age)[1], max = HDInterval::hdi(age)[2])

  stan_data <- posterior::as_draws_df(stan_output) %>%
    select(starts_with('modelled_age'))

  names(stan_data) <- x$label

  age_distrib <- stan_data %>% pivot_longer(everything(), names_to = 'label', values_to = 'age') %>%
    inner_join(select(x, label, depth), by = 'label')

  output <- list(
    input = x,
    stan_output = stan_output,
    modelled_ages = output_ages,
    modelled_ages_summary = output_ages_summary,
    output_ages = output_ages_summary,
    age_distrib = age_distrib,
    model_family = 'uniform',
    implementation_family = implementation_family,
    implementation = paste0('uniform_', implementation_family)
  )

  class(output) <- append(class(output), 'agedepthmodel')

  output
}


run_uniform_simple <- function(x, chains = 4, warmup = 1500, iter = 4000) {

  stan_model_data <- list(
    N = nrow(x),
    depth = x$depth,
    age = x$calib_date_mean,
    age_sd = x$calib_date_sd
  )

  model <- stanmodels$uniformsimple

  stan_output <- rstan::sampling(
    object = model,
    data = stan_model_data,
    chains = chains,
    warmup = warmup,
    iter = iter,
    cores = parallel::detectCores()
  )

  uniform_postprocess(x, stan_output, 'simple')

  }

run_uniform_full <- function(x, chains = 4, warmup = 1000, iter = 3000) {
  distribs <- x$calibrated_dist
  num_dates <- nrow(x)
  num_probs <- max(vapply(distribs, nrow, numeric(1)))

  ages_matrix <- matrix(nrow = num_probs, ncol = num_dates)
  probs_matrix <- matrix(nrow = num_probs, ncol = num_dates)

  for(i in 1:num_dates) {
    dates<-distribs[[i]]$cal_bp
    probs<-distribs[[i]]$prob
    shortfall <- num_probs - length(probs)
    ages_matrix[,i] = c(rep(-1,shortfall), dates)
    probs_matrix[,i] = c(rep(0,shortfall), probs)
  }

  cc <- as.data.frame(ccurve(resample=10))
  names(cc) <- c('V1', 'V2', 'V3')



  stan_model_data <- list(
    num_dates = num_dates,
    num_probs = num_probs,
    ages = t(ages_matrix),
    probs = t(probs_matrix),
    depth = x$depth,
    calcurve_size = nrow(cc),
    calibrated_ages = cc$V1,
    c14_ages = cc$V2,
    cc_sd = cc$V3
  )

  init_values <- function() {
    list(det_ages_dist = x$calib_date_mean)
  }

  model <- stanmodels$uniformactualdistributions

  stan_output <- rstan::sampling(
    object = model,
    chains = chains,
    data = stan_model_data,
    warmup = warmup,
    iter = iter,
    cores = parallel::detectCores(),
    init = init_values
  )

  uniform_postprocess(x, stan_output, 'full')
}

run_sequential_full <- function(x, chains = 4, warmup = 1000, iter = 3000) {
  distribs <- x$calibrated_dist
  num_dates <- nrow(x)
  num_probs <- max(vapply(distribs, nrow, numeric(1)))

  ages_matrix <- matrix(nrow = num_probs, ncol = num_dates)
  probs_matrix <- matrix(nrow = num_probs, ncol = num_dates)

  for(i in 1:num_dates) {
    dates<-distribs[[i]]$cal_bp
    probs<-distribs[[i]]$prob
    shortfall <- num_probs - length(probs)
    ages_matrix[,i] = c(rep(-1,shortfall), dates)
    probs_matrix[,i] = c(rep(0,shortfall), probs)
  }


  stan_model_data <- list(
    num_dates = num_dates,
    num_probs = num_probs,
    ages = t(ages_matrix),
    probs = t(probs_matrix)
  )

  model <- stanmodels$sequentialactualdistributions

  stan_output <- rstan::sampling(
    object = model,
    chains = chains,
    data = stan_model_data,
    warmup = warmup,
    iter = iter,
    cores = parallel::detectCores()
  )

  stan_data <- posterior::as_draws_df(stan_output) %>%
    select(starts_with('true_ages'))

  names(stan_data) <- x$label

  age_distrib <- stan_data %>% pivot_longer(everything(), names_to = 'label', values_to = 'age') %>%
    inner_join(select(x, label, depth), by = 'label')

  output <- list(input = x,
    stan_output = stan_output,
    age_distrib = age_distrib)

  class(output) <- append(class(output), 'agedepthmodel')

  output
}

run_sequential_simple <- function(x, chains = 4, warmup = 1000, iter = 3000) {

  stan_model_data <- list(
    N = nrow(x),
    dates = x$calib_date_mean,
    error = x$calib_date_sd
  )

  model <- stanmodels$sequentialsimple

  stan_output <- rstan::sampling(
    object = model,
    data = stan_model_data,
    warmup = warmup,
    iter = iter,
    chains = chains,
    cores = parallel::detectCores()
  )

  stan_data <- posterior::as_draws_df(stan_output) %>%
    select(starts_with('true_ages'))

  names(stan_data) <- x$label

  age_distrib <- stan_data %>% pivot_longer(everything(), names_to = 'label', values_to = 'age') %>%
    inner_join(select(x, label, depth), by = 'label')

  output <- list(input = x,
    stan_output = stan_output,
    age_distrib = age_distrib)

  class(output) <- append(class(output), 'agedepthmodel')

  output
}


#' @title Run an age-depth model
#'
#' @description Runs the specified age-depth model and returns an object of type `agedepthmodel`
#'
#' The input data frame should contain four columns: `label`, a label for the observation, `age`, an uncalibrated C14 age,
#' `error`, the error of the C14 age and `depth`, the depth at which the observation was taken.
#'
#' The model parameter can take one of the values bacon_simple, uniform_simple, sequential_simple, sequential_full, uniform_full.
#'
#' @param x  A data frame containing age depth data
#' @param model A string specifying the model to be run
#' @param chains The number of MCMC chains to be run
#' @param warmup The number of MCMC iterations to be discarded as warmup
#' @param iter The total number of MCMC iterations
#' @param modelled_points (Bacon only) the number of points to model. The number of calculated rates is one less that this number.
#' @param output_points (Bacon only) the number of depths to provide an age for in the output
#' @param a_alpha (Bacon only) the a input for the alpha (rate) prior
#' @param b_alpha (Bacon only) the b input for the alpha (rate) prior
#' @param a_w (Bacon only) the a input for the memory prior
#' @param b_w (Bacon only) the b input for the memory prior
#' @param t_df (Bacon only) the number of degrees of freedom for the t-distribution in the likelihood function
#' @param use_normal (Bacon only) If 1, use the normal distribution instead of the t-distribution in the likelihood
#'
#' @examples
#' # run_model(sample_data, 'bacon_simple', modelled_points = 7, a_alpha = 140, b_alpha = 0.1, a_w = 5, b_w = 5, t_df = 3)
#' # run_model(MSB2K, 'bacon_simple', modelled_points = 21, a_alpha = 1.5, b_alpha = 0.075, a_2 = 5, b_w = 5, t_df = 3)
#'
#' @export
run_model <- function(x, model = c('uniform_simple', 'bacon_simple', 'sequential_simple', 'sequential_full', 'uniform_full', 'bacon_full'), ...) {
  model <- match.arg(model)
  x <- preprocess_input(x)
  if (model == 'bacon_simple') {
    run_bacon_simple(x, ...)
  } else if (model == 'uniform_simple') {
    run_uniform_simple(x, ...)
  } else if (model == 'sequential_simple') {
    run_sequential_simple(x, ...)
  } else if (model == 'bacon_full') {
    stop('bacon_full model not supported in this package.')
    run_bacon(x, ...)
  } else if (model == 'sequential_full') {
    run_sequential_full(x, ...)
  } else if (model =='uniform_full') {
    run_uniform_full(x, ...)
  }

}
