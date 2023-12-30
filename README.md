# Bayesian age depth models in stan

## Background

These files implement age-depth models in stan, prepared as part of my dissertation for the University of Sheffield's MSc in Statistics.

 - OxCal's sequential model (see https://c14.arch.ox.ac.uk/oxcal.html)
 - OxCal's uniform deposition model
 - Bacon (see https://github.com/cran/rbacon)
 

In the case of the Oxcal models, a version using both a normal approximation to calibrated C14 dates, as well as a full version of the C14
calibrated date distribution. The operable bacon model only implements the normal approximation. A full version of my implementation of the 
bacon model is in the extras/ directory, but that requires the CmdStanR package. The R code needed to run the full bacon version is in the 
run_model.R file, but the run_model() function has a call to stop() to prevent it from running. Nevertheless, it is possible to inspect the 
stan code.


## Arrangement of code

R code to call the stan functions is in R/run_model.R. R code to use the objects created by run_model to predict ages and to plot age-depth
models is in R/tools.R

Stan code for the majority of models is in inst/stan.

The code for the implementation of bacon using full C14 distributions is in extras/. That directory also contains code to run the Oxcal 
uniform model. To run the OxCal uniform model, OxCal must be installed and the call to setOxcalExecutablePath updated appropriately.

## Usage

In principle it should be possible to install this as a package by first installing the devtools package and then running devtools::install_github('https://github.com/stanagedepth/stanagedepth'). This can be a very slow process (because of the 
need to recompile the stan files), and sometimes somewhat unreliable. However if it works you will end up with an 
installed R library called stanagedepth with documented functions:

 - sample_data() will return a dataframe containing the sample data used in my "usage" chapter.
 - run_model() takes a data frame and a string specifying the model to be run: one of sequential_simple, uniform_simple, bacon_simple, sequential_full, uniform_full or (not working directly) bacon_full.

Alternatively, it should be possible to install the stan and R files separately. It would be necessary to replace the references to the 
stanmodels object with a call to rstan::stan_model to compile and load the models.
