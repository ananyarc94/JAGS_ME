rm(list = ls())
set.seed(100)

# Load packages
#library(mgcv)
library(ggplot2)
library(gridExtra)
#library(lme4)
library(rjags) # must have JAGS installed on computer to use this package
library(coda)
library(mvtnorm)
library(tidymodels)
library(dplyr)

# Source functions
function_directory <- getwd()
source(paste0(function_directory, "/code/process_data.R"))
source(paste0(function_directory, "/code/plot_gam.R"))

#################################################################################
# Define variable and process data for analysis
#################################################################################

#### Old data analysis 8/25/21
# Load the data
sim_data <- read.csv(paste0(function_directory, "/code/simulated_data_2021.08.25.csv"))

# Fix the fatigue values that differ from updated dataset
sim_data.new <- read.csv(paste0(function_directory, "/code/simulated_data_2021.12.21.csv"))
sim_data$fatigue_mean_intensity_score <- sim_data.new$fatigue_mean_intensity_score
sim_data$fatigue_mean_interference_score <- sim_data.new$fatigue_mean_interference_score

# Create variable for mean activity bout duration
## activity threshold = 10 counts/min
sim_data$mean_activity_bout_duration10 <- 1 / sim_data$astp_as_is_10_orig

# Discard missing response values
sim_data <- sim_data[which(!is.na(sim_data[ , c("fatigue_mean_intensity_score")])),]

# Vector of variables to take difference from baseline
vars <- c("fatigue_mean_intensity_score", "fatigue_mean_interference_score",
          "mvpa_avg_as_is_1952_orig", "mean_activity_bout_duration10")

# Variable indicating group assignment
group <- "Randomization"

# Pre-process the data
data <- preprocess_data(sim_data, vars, group)


# Vectors of all the covariates for analysis
## response
y <- "fatigue_mean_intensity_score_diff"
y_sc <- "fatigue_mean_intensity_score_diff_sc"
# y <- "fatigue_mean_interference_score_diff"
# y_sc <- "fatigue_mean_interference_score_diff_sc"

## physical activity vars
pa_vars <- c("mvpa_avg_as_is_1952_orig_diff_sc", "mean_activity_bout_duration10_diff_sc")
plot_title <- c("MVPA", "Mean ABD") # titles for plotting

## non-PA covariates
covars <- c("age100", "BMI_obese", "ht_sample", "breast_sample")
covars <- c(paste0(covars, ".C"), paste0(covars, ".I"))



#################################################################################
# Build the model
#################################################################################
#convert the categorical covariates into dummy variables

rec_obj <- recipe(fatigue_mean_intensity_score_diff_sc ~ ., data = data)
ind_vars <- rec_obj %>%
  step_dummy(all_nominal_predictors())
trained_rec <- prep(ind_vars, training = data)
train_data <- bake(trained_rec, new_data = data)

#selecting the relevant columns to form the design matrix X
X.mat <- train_data[, c(21,23,24,34:50)]



model_string <- "model{
  for (i in 1:n){

    # response model
    mu[i] <- beta%*%X.mat[i,] + g*z[i]
    y[i] ~ dnorm (mu[i], tau.y)


    # MVPA model
    z[i] ~ dnorm(0, taue)
    # classical measurement error model
    mvpa_avg_as_is_1952_orig_diff_sc[i] ~ dnorm(z[i], tauu)



  }

  # prior distributions
  beta ~ dmnorm(rep(0, p), prec)
  g ~  dnorm(0, 0.1)




  # prior distributions error model
  tau.y ~ dgamma(0.01, 0.01)
  tauu ~ dgamma(0.01, 0.01)
  taue ~  dgamma(0.01, 0.01)




}"



## Build the dataset to fit JAGS model
data.jags <- list(y = as.vector(data$fatigue_mean_intensity_score_diff_sc), n = length(data$fatigue_mean_intensity_score_diff_sc),
                  mvpa_avg_as_is_1952_orig_diff_sc = as.vector(data$mvpa_avg_as_is_1952_orig_diff_sc),
                  X.mat = X.mat, p = ncol(X.mat), prec = diag(rep(0.01,ncol(X.mat))))

# Intilialize JAGS model
jags.mod <- jags.model(textConnection(model_string), data = data.jags,
                       inits = list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 2021), # this sets the seed in JAGS to get the same results
                       n.chains = 1, n.adapt = 1000)

# Specify number of burn-in
update(jags.mod, n.burn = 5000)

# Obtain JAGS samples
samples.jags <- jags.samples(jags.mod, c("beta", "g","tauu", "tau.y"), n.iter = 50000, thin = 50)
samples.jags

## Obtain posterior samples as a mcmc.list object
samples.coda <- coda.samples(jags.mod, c("beta", "g","tauu", "tau.y"), n.iter = 50000, thin = 50)
samples.array <- as.array(samples.coda)
dimnames(samples.array)[[2]] <- c(colnames(X.mat),"g","tauu","tau.y")


## Traceplots for the parameter estimates
par(mfrow = c(3, 3))
for(i in 1:ncol(samples.array)){
  traceplot(mcmc(samples.array[,i]),main = dimnames(samples.array)$var[i])
}
#traceplot(mcmc(samples.array[,ncol(samples.array)]), main = "tauu")

#ACF plot for gamma
acf(samples.array[,(ncol(samples.array)-1)], main = "g")
acf(samples.array[,(ncol(samples.array))], main = "tauu")

## effect sample sizes
par(mfrow=c(1,1))
Neff <- effectiveSize(samples.coda)
plot(1:length(Neff), Neff, main = "Effective sample size", xlab = "parameter")



beta.list <- lapply(1:20, function(i) {
  ggplot(data = as.data.frame(samples.array[ ,i]), aes(x = 1:1000, y = samples.array[ ,i] ))+
    geom_line() + xlab("Index") +  theme_bw()+
    ggtitle(colnames(samples.array)[i])+ylab(NULL)
})

png(file="./plots/1var/beta_1.png")
ggarrange(plotlist = beta.list, nrow = 4, ncol = 3)$`1`
dev.off()

png(file="./plots/1var/beta_2.png")
ggarrange(plotlist = beta.list, nrow = 4, ncol = 3)$`2`
dev.off()

g <-  ggplot(data = as.data.frame(samples.array[ ,21]), aes(x = 1:1000, y = samples.array[ ,21] ))+
    geom_line() + xlab("Index")+  theme_bw() +
    ggtitle(colnames(samples.array)[21])+ylab(NULL)

tau.u <-  ggplot(data = as.data.frame(samples.array[ ,22]), aes(x = 1:1000, y = samples.array[ ,22] ))+
  geom_line() + xlab("Index")+  theme_bw() + ylim(0,50)+
  ggtitle(colnames(samples.array)[22])+ylab(NULL)

tau.y <- ggplot(data = as.data.frame(samples.array[ ,23]), aes(x = 1:1000, y = samples.array[ ,23] ))+
  geom_line() + xlab("Index")+  theme_bw() + ylim(0,10)+
  ggtitle(colnames(samples.array)[23])+ylab(NULL)

png(file="./plots/1var/g_and_sigma_plots.png")
grid.arrange(g, tau.u, tau.y, nrow = 2, ncol = 2)
dev.off()


