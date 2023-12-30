library(dplyr)
library(tidyr)
library(posterior)
library(stringr)
library(bayesplot)
library(oxcAAR)
library(V8)


setOxcalExecutablePath('/cloud/project/OxCal/bin/OxCalLinux')

ox_data <- read.csv('/cloud/project/oxcal_data.csv')
ox_data <- preprocess_input(ox_data)

unif_oxcal <- function(x) {
  # See https://c14.arch.ox.ac.uk/oxcalhelp/hlp_analysis_oper.html



  oxcal_script <- paste0(
    'Options() {Ensembles = 50;};
                 U_Sequence("U_SEQUENCE",5) {
                      Boundary("Bottom_boundary"){};',
    paste0((x %>%
              arrange(desc(depth)) %>%
              mutate(row_num = seq(1:nrow(.)), r_date = paste0('R_date("', label, '", ', age, ', ', error, ') { z=', depth, ';};')))$r_date, collapse = ''),
    'Boundary("Top_boundary") {};
                      };'
  )

  oxcal_output_file <- executeOxcalScript(oxcal_script)
  oxcal_output_text <- readOxcalOutput(oxcal_output_file)

  #The standard library's oxcal parser won't help. Need to use a javascript library to execute the OxCal output
  js_context <- v8()
  # We need to create the objects in our javascript context that the Oxcal output data will fill up
  js_context$eval("ocd={}; calib={}; model={};")
  js_context$eval(oxcal_output_text)
  oxcal_ocd <- js_context$get('ocd')
  oxcal_calib <- js_context$get('calib')
  oxcal_model <- js_context$get('model')

  #Age values from the model ensemble
  age_values_matrix <- data.frame(filter(oxcal_model$element, age_depth_ensembles != 'NULL')$age_depth_ensembles[[1]])

  names(age_values_matrix) <- filter(oxcal_model$element, age_depth_z != 'NULL')$age_depth_z[[1]]
  age_values <- age_values_matrix %>%
    pivot_longer(everything(), names_to = 'depth', values_to = 'age') %>%
    mutate(depth = as.numeric(depth)) %>%
    inner_join(select(x, label, depth), by = 'depth')  %>%
    mutate(age = 1950 - age)

  rates <- oxcal_model$element %>%
    filter(age_depth_rate != 'NULL') %>%
    select(c('age_depth_rate', 'age_depth_ratesd')) %>%
    mutate(age_depth_rate = unlist(age_depth_rate)[1], age_depth_ratesd = unlist(age_depth_ratesd)[1])

  N <- 500 # number of outputs

  age_depth <- data.frame(list(
    iter = numeric(0),
    depth = numeric(0),
    age = numeric(0)
  ))

  ratesamples <- rnorm(N, rates$age_depth_rate[1], rates$age_depth_ratesd[1])
  for (i in 1:N) {
    rate = ratesamples[i]
    base_age = sample((age_values %>% filter(depth == x$depth[1]))$age, size = 1)
    base_depth = x$depth[1]
    o <- data.frame(iter = rep(i, nrow(x)), depth = x$depth, age = numeric(nrow(x)))
    o <- o %>%
      mutate(age = base_age + (depth - base_depth) * rate)
    age_depth <- rbind(age_depth, o)
  }


  age_depth_summary <- age_depth %>%
    select(-iter) %>%
    group_by(depth) %>%
    summarise(
      mean = mean(age),
      min = HDInterval::hdi(age)[1],
      max = HDInterval::hdi(age)[2]
    )

  return(list(
    ages_values = age_values,
    rates = rates,
    age_depth = age_depth,
    age_depth_summary = age_depth_summary
  ))
}
