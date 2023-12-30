library(ggplot2)
library(ggridges)


#' @title Predict age for a given depth for an age-depth model
#'
#' @description Use the input value to determine the age at a given depth.
#' @param x An object of type `agedepthmodel`
#' @param depth The depth at which the age is required
#' @return A named vector containing the mean age according to the model, together with the minimum and maximum values (95% HDI)
#' @export
predict.agedepthmodel <- function(x, depth) {
  if(!'agedepthmodel' %in% class(x))
    stop('Invalid argument class')
  if(is.null(x$output_ages)) {
    warning('predict.agedepthmodel: cannot predict an age if no age-depth model is produced. Try, e.g. uniform deposition or bacon model')
    return(NULL)
  }
  if(depth < min(x$output_ages$depth) || depth > max(x$output_ages$depth)) {
    warning(paste0('predict.agedepthmodel: depth ', depth, ' out of model range'))
  }

  min  <- approx(x$output_ages$depth, x$output_ages$min, depth)$y
  mean <- approx(x$output_ages$depth, x$output_ages$mean, depth)$y
  max  <- approx(x$output_ages$depth, x$output_ages$max, depth)$y

  c(min = min, mean = mean, max = max)
}

#' @title Plot one or two age-depth models
#'
#' @description Plot of one or two age-depth models, including input and posterior ages where available.
#' @param x An object of type `agedepthmodel`
#' @param y A second, optional, object of type `agedepthmodel`
#' @param colours The colours to be used for the two age-depth models
#' @param ageunits The description of the units to be used on the age axis
#' @param depthunits The description of the units to be used on the depth axis
#' @param include_input_ages Boolean. If true, plots the ages used as input to the model
#' @param include_labels Boolean. If true, include the observation labels from the input data
#' @return An object of type `ggplot`
#' @examples
#' # plot(sample_data())
#'
#' @export
plot.agedepthmodel <- function(x,
                               y = NULL,
                               colours = c('red', 'green'),
                               ageunits = 'years bp',
                               depthunits = 'cm',
                               include_input_ages = TRUE,
                               include_labels = TRUE) {
  if(!('agedepthmodel' %in% class(x) || (!is.null(y) && !'agedepthmodel' %in% class(y))))
    stop('Invalid argument classes')

  if(is.null(x$input))
    stop('First argument must include the input')

  ages_original <- x$input %>% select(label, depth, calibrated_dist) %>% unnest(cols = c(calibrated_dist))

  include_x_agedepth <- !is.null(x$output_ages)
  include_y_agedepth <- !is.null(y) && !is.null(y$output_ages)
  include_x_agedistrib <- !is.null(x$age_distrib)
  include_y_agedistrib <- !is.null(y) && !is.null(y$age_distrib)

  if(include_labels) {
    labels <- select(ages_original, c('label', 'depth', age = 'cal_bp'))

    if(include_x_agedistrib) {
      labels <- labels %>%
        bind_rows(select(x$age_distrib, c('label', 'depth', 'age')))
    }

    if(include_y_agedistrib) {
      labels <- labels %>%
        bind_rows(select(y$age_distrib, c('label', 'depth', 'age')))
    }


    labels <- labels %>%
      group_by(label, depth) %>%
      summarise(age = min(age), .groups = 'keep')
  }

    agedepth <- ggplot() +
      {if(include_x_agedepth) geom_line(data = x$output_ages, aes(y = depth, x = mean), colour = colours[1], alpha = 0.75)} +
      {if(include_x_agedepth) geom_line(data = x$output_ages, aes(y = depth, x = max), colour = colours[1], linetype = 'dashed', alpha = 0.5)} +
      {if(include_x_agedepth) geom_line(data = x$output_ages, aes(y = depth, x = min), colour = colours[1], linetype = 'dashed', alpha = 0.5)} +

      {if(include_y_agedepth) geom_line(data = y$output_ages, aes(y = depth, x = mean), colour = colours[2], alpha = 0.75)} +
      {if(include_y_agedepth) geom_line(data = y$output_ages, aes(y = depth, x = max), colour = colours[2], linetype = 'dashed', alpha = 0.5)} +
      {if(include_y_agedepth) geom_line(data = y$output_ages, aes(y = depth, x = min), colour = colours[2], linetype = 'dashed', alpha = 0.5)} +

      {if(include_x_agedistrib) geom_density_ridges(colour = colours[1], fill = colours[1], alpha = 0.1, rel_min_height = 0.001, data = x$age_distrib, aes(y = depth, x = age, group = label))} +

      {if(include_y_agedistrib) geom_density_ridges(colour = colours[2], fill = colours[2], alpha = 0.1, rel_min_height = 0.001, data = y$age_distrib, aes(y = depth, x = age, group = label))} +

      {if(include_input_ages) geom_ridgeline(colour = 'grey', fill = 'grey', alpha = 0.1, data = ages_original, aes(x = cal_bp, y = depth, height = 3 * prob, group = label))} +

      {if(include_labels) geom_text(data = labels, mapping = aes(x = age, y = depth + 0.005 - 0.005 - 0.005, label = label), hjust = 0)} +

      theme_bw() +
      xlab(paste0('Age / ', ageunits)) +
      ylab(paste0('Depth / ', depthunits)) +
      scale_y_reverse()

  agedepth

}
