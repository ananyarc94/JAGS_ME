################################################################################
# This program is an implementation of a measurement error model either linear or
# splines using JAGS.
################################################################################


rm(list = ls())
set.seed(100)

# Load necessary packages
library(mgcv)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(rjags)
library(coda)
library(ggpubr)
#library(bSpline)
# Source functions
function_directory <- getwd()
source(paste0(function_directory, "/code/process_data.R"))

#################################################################################
# Pre-process the simulated data
#################################################################################

# Load the data
sim_data <- read.csv(paste0(function_directory, "/code/simulated_data_2021.08.25.csv"))

#################################################################################
# Implement the linear model
#################################################################################

# Preprocess the data
data_processed = dataPreprocess(sim_data)

X.mat = data_processed$X.mat
MVPA_C = data_processed$MVPA_C
MVPA_I = data_processed$MVPA_I
ABD_C = data_processed$ABD_C
ABD_I = data_processed$ABD_I

## Build the dataset to fit JAGS model
data.jags <- list(y = as.vector(data$fatigue_mean_intensity_score_diff_sc), n = length(data$fatigue_mean_intensity_score_diff_sc),
                  MVPA_C = as.vector(MVPA_C),MVPA_I = as.vector(MVPA_I), ABD_C = as.vector(ABD_C), ABD_I = as.vector(ABD_I),
                  X.mat = X.mat, p = ncol(X.mat), prec = diag(rep(0.01,ncol(X.mat))))

# Intilialize JAGS model
jags.mod <- jags.model(textConnection(model_string), data = data.jags,
                       inits = list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 2021), # this sets the seed in JAGS to get the same results
                       n.chains = 1, n.adapt = 1000)

# Specify number of burn-in
update(jags.mod, n.burn = 5000)

# Obtain JAGS samples
samples.jags <- jags.samples(jags.mod, c("beta", "g1", "g2", "g3","g4","sigma_u1","sigma_u2",
