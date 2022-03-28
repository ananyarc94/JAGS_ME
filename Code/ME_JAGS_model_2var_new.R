################################################################################
# This program creates a factor-by-curve interaction model on the simulated
# data by creating a Bayesian graphical model with JAGS.
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
library(recipes)
library(ggpubr)
# Source functions
function_directory <- setwd("C:/Users/anany/Desktop/JAGSME")
source(paste0(function_directory, "/code/process_data.R"))

#################################################################################
# Pre-process the simulated data
#################################################################################

# Load the data
sim_data <- read.csv(paste0(function_directory, "/code/simulated_data_2021.08.25.csv"))

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
pa_vars_orig <- c("mvpa_avg_as_is_1952_orig_diff", "mean_activity_bout_duration10_diff")
plot.title <- c("MVPA", "Mean ABD") # titles for plotting

## non-PA covariates
covars <- c("age100", "BMI_obese", "ht_sample", "breast_sample")
covars <- c(paste0(covars, ".C"), paste0(covars, ".I"))

#################################################################################
# Change order of data to be blocks for each group and time period. This makes
# plotting curve estimated in JAGS easier. New data will be control group at M3,
# control group at M6, ..., treatment group at month 6, treatment group at
# month 12.
################################################################################

## Separate data by group (control vs intervention)
indC <- which(data$Randomization == 0)
indI <- which(data$Randomization == 1)
dataC <- data[indC, ]
dataI <- data[indI, ]

### Get indices for time periods for each group
# Month 3
t1.C <- which(data$svy_interval[indC] == 1)
t1.I <- which(data$svy_interval[indI] == 1)

# Month 6
t2.C <- which(data$svy_interval[indC] == 2)
t2.I <- which(data$svy_interval[indI] == 2)

# Month 12
t3.C <- which(data$svy_interval[indC] == 3)
t3.I <- which(data$svy_interval[indI] == 3)

### Reorder response
yC <- as.matrix(dataC[ , y_sc], ncol = 1)
yI <- as.matrix(dataI[ , y_sc], ncol = 1)
colnames(yC) <- colnames(yI) <- y_sc
Y <- as.matrix(c(yC[c(t1.C, t2.C, t3.C)], yI[c(t1.I, t2.I, t3.I)]), ncol = 1)

### Reorder matrix for PA variables
xC <- as.matrix(dataC[ , pa_vars], ncol = length(pa_vars))
xI <- as.matrix(dataI[ , pa_vars], ncol = length(pa_vars))
X <- rbind(xC[c(t1.C, t2.C, t3.C),], xI[c(t1.I, t2.I, t3.I),])

### Reorder matrix for other covariates
zC <- as.matrix(dataC[ , covars[grepl(".C", covars, fixed = TRUE)]]) # made change here!!!!
zI <- as.matrix(dataI[ , covars[grepl(".I", covars, fixed = TRUE)]])
Z <- cbind(rbind(zC[c(t1.C, t2.C, t3.C),], matrix(0, nrow = nrow(zI), ncol = ncol(zC))),
           rbind(matrix(0, nrow = nrow(zC), ncol = ncol(zI)), zI[c(t1.I, t2.I, t3.I),]))

## Build a dataset containing variables of interest
int.and.group <- rbind(dataC[c(t1.C, t2.C, t3.C), c("int.time", group)],
                       dataI[c(t1.I, t2.I, t3.I), c("int.time", group)])
data.sim <- data.frame(Y, Z, int.and.group, X)
colnames(data.sim) <- c(y_sc, colnames(zC), colnames(zI), "int.time", group, pa_vars)

# Change age variables to numeric
# (turns into factor when combined into matrix with other factor variables)
data.sim[, covars[grepl("age", covars)]] <- apply(data.sim[, covars[grepl("age", covars)]], 2, as.numeric)
#convert the categorical covariates into dummy variables



#################################################################################
# Build the model
#################################################################################
# #convert the categorical covariates into dummy variables
rec_obj <- recipe(fatigue_mean_intensity_score_diff_sc ~ ., data = data.sim) # made change here!!!
ind_vars <- rec_obj %>%
  step_dummy(all_nominal_predictors())
trained_rec <- prep(ind_vars, training = data.sim)
train_data <- bake(trained_rec, new_data = data.sim)

# #Creating the column for interval 1
# int2 = which(X.mat$int.time_X2 == 1)
# int3 = which(X.mat$int.time_X3 == 1)
# int4 = which(X.mat$int.time_X4 == 1)
# int5 = which(X.mat$int.time_X5 == 1)
# int6 = which(X.mat$int.time_X6 == 1)
# int1 = setdiff(1:275, sort(c(int2, int3, int4, int5, int6)))

int.time_X1 = rep(0, nrow(train_data))
int.time_X1[int1] = 1

#selecting the relevant columns to form the design matrix X
# X.mat <- data.frame(int.time_X1,train_data[, c(34:38, 23,24, 39:50)])
X.mat <- data.frame(int.time_X1,train_data[, c(19:23, 1, 2, 7:18)])

#splitting up the MVPA and ABD into the control and treatment groups
MVPA_C <- ifelse(train_data$Randomization == 0, train_data$mvpa_avg_as_is_1952_orig_diff_sc, 0)
MVPA_I <- ifelse(train_data$Randomization == 1, train_data$mvpa_avg_as_is_1952_orig_diff_sc, 0)
ABD_C <- ifelse(train_data$Randomization ==0, train_data$mean_activity_bout_duration10_diff_sc,0)
ABD_I <- ifelse(train_data$Randomization ==1, train_data$mean_activity_bout_duration10_diff_sc,0)


#################################################################################
# Define JAGS model
#################################################################################

model_string <- "model{
  for (i in 1:n){

    # response model
    mu[i] <- beta%*%X.mat[i,] + g1*z1[i] + g2*z2[i] +g3*z3[i] +g4*z4[i]
    y[i] ~ dnorm (mu[i], tau.y)


    # MVPA model for control group
    z1[i] ~ dnorm(0, taue1)
    # classical measurement error model
    MVPA_C[i] ~ dnorm(z1[i], tauu1)

    # MVPA model for treatment group
    z2[i] ~ dnorm(0, taue2)
    # classical measurement error model
    MVPA_I[i] ~ dnorm(z2[i], tauu2)

    # MABD model for control group
    z3[i] ~ dnorm(0, taue3)
    # classical measurement error model
    ABD_C[i] ~ dnorm(z3[i], tauu3)

    # MABD model for treatment group
    z4[i] ~ dnorm(0, taue4)
    # classical measurement error model
    ABD_I[i] ~ dnorm(z4[i], tauu4)

  }

  # prior distributions
  beta ~ dmnorm(rep(0, p), prec)
  g1 ~  dnorm(0, 0.1)
  g2 ~  dnorm(0, 0.1)
  g3 ~  dnorm(0, 0.1)
  g4 ~  dnorm(0, 0.1)



  # prior distributions error model
  tau.y ~ dgamma(0.01, 0.01)
  tauu1 ~ dgamma(0.01, 0.01)
  tauu2 ~ dgamma(0.01, 0.01)
  tauu3 ~ dgamma(0.01, 0.01)
  tauu4 ~ dgamma(0.01, 0.01)
  taue1 ~  dgamma(0.01, 0.01)
  taue2 ~  dgamma(0.01, 0.01)
  taue3 ~  dgamma(0.01, 0.01)
  taue4 ~  dgamma(0.01, 0.01)

  # converting the precisions back to the variances
  sigma_y <- 1/tau.y
  sigma_u1 <- 1/tauu1
  sigma_u2 <- 1/tauu2
  sigma_u3 <- 1/tauu3
  sigma_u4 <- 1/tauu4

}"

#################################################################################
# Build the JAGS model
#################################################################################

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
                                         "sigma_u3","sigma_u4","sigma_y"), n.iter = 50000, thin = 50)
samples.jags

################################################################################
# Organize the results
################################################################################

samples <- c(apply(samples.jags$b, 1, mean),apply(samples.jags$g1, 1, mean),apply(samples.jags$g2, 1, mean),
             apply(samples.jags$g3, 1, mean),apply(samples.jags$g4, 1, mean))
W_mat <- cbind(MVPA_C, MVPA_I, ABD_C, ABD_I)

# Obtain index of covariates and intercept
cov <- 1:ncol(X.mat)

# Get mean and sd from MCMC samples
theta.sample <- samples.jags$beta
est <- apply(theta.sample, 1, mean)
se <- apply(theta.sample, 1, sd)

# Calculate p-value for whether estimate != 0
pvals <-  2 * (1 - pnorm(abs(est / se)))

## Combine output and print results
results <- cbind(est, se, pvals)
rownames(results) <- colnames(X.mat)
print(results, digits = 2)


## Organize the curve output

# Combine scaled PA vars into matrix
X <- data.sim[ , pa_vars]

# Combine original units of PA vars into matrix
X.orig <- rbind(dataC[c(t1.C, t2.C, t3.C), pa_vars_orig],
                dataI[c(t1.I, t2.I, t3.I), pa_vars_orig])
# List to save curve results
curve <- list()

# Loop over each smooth function
# (one curve for each group and for each PA variable)
for(i in 1:(2*length(pa_vars))){

  # Data is grouped as MVPA control group, MVPA trt group, ABD control group, ABD trt group
  # Specify type = 0 for control in data groups 1, 3
  if(i %in% c(1, 3)){
    type = 0
  }else{
    # Specify type = 1 for treatment in data groups 2, 4
    type = 1
  }
  # Indices of corresponding ME variables
  ind.spl <- ncol(X.mat) + i


  # Calculate estimated curve for each MCMC sample
  #fHat.all <- rep(apply(X.mat,2,mean)%*% apply(samples.jags$beta, 1, mean),length(MVPA_C)) + MVPA_C * rep(apply(samples.jags$g1, 1, mean),length(MVPA_C))
  #          + rep(mean(MVPA_I)*apply(samples.jags$g2, 1, mean),length(MVPA_C)) + rep(mean(ABD_C)*apply(samples.jags$g3, 1, mean),length(MVPA_C)) + rep(mean(ABD_I)*apply(samples.jags$g4, 1, mean),length(MVPA_C))

  fHat.all <- W_mat[,i] * samples[ind.spl]

  # Backtransform estimated curves into original units
  fHat.tilde <- fHat.all * sd(data[ , y]) + mean(data[ , y])

  # Obtain mean and 95% CI for curve
  fHat.mean <- mean(fHat.tilde)
  #lower <- quantile(fHat.tilde, 0.025)
  #upper <- quantile(fHat.tilde, 0.975)


  # Obtain x values for MVPA for plotting
  # data groups 1, 2 correspond to MVPA
  if(i %in% c(1, 2)){
    x.orig <- X.orig[, 1]
    x <- X[ , 1]
    pa <- "MVPA"
  }else{     # Obtain x values for ABD for plotting
    # data groups 3, 4 correspond to ABD
    x.orig <- X.orig[, 2]
    x <- X[ , 2]
    pa <- "ABD"
  }

  # Combine x and y values into list
  curve[[paste0(pa, i)]] <- data.frame(x.orig = x.orig, s.x = fHat.tilde)

  # Smooth CI estimates and save into list
  #curve[[paste0(pa, i)]]$lower <- lm(lower ~ x)$fitted.values
  #curve[[paste0(pa, i)]]$upper <- lm(upper ~ x)$fitted.values

}

# Obtain intercepts in original scale
ni <- length(unique(data$int.time))
int <- est[1:ni] * sd(data[ , y]) + mean(data[ , y])

# Obtain curves with intercepts for each group
curve.with.int <- list()

for(pa in names(curve)){
  # Add time-dependent intercept to smooth functions for each group
  ## If in control group (data groups 1, 3)
  if(any(grepl('1', strsplit(pa, '')[[1]])) | any(grepl('3', strsplit(pa, '')[[1]]))){
    # Add time-and-group-dependent interept to each curve (mean and CI)
    # Save final curve with intercept into dataframe
    curve.with.int[[pa]] <- lapply(1:(length(int)/2), function(ind){x = int[ind];
    int.ind <- 1:nrow(dataC);
    data.frame(x = curve[[pa]]$x.orig[int.ind],
               s.x = x + curve[[pa]]$s.x[int.ind],
               #lower = x + curve[[pa]]$lower[int.ind],
               #upper = x + curve[[pa]]$upper[int.ind],
               type = ind);})
    # Label as control group
    ctrl <- TRUE
  }
  ## If in treatment group (data groups 2, 4)
  else{
    # Add time-and-group-dependent interept to each curve (mean and CI)
    # Save final curve with intercept into dataframe
    curve.with.int[[pa]] <- lapply((length(int)/2 +1):length(int), function(ind){x = int[ind];
    int.ind <- (nrow(dataC)+1):nrow(data.sim);
    data.frame(x = curve[[pa]]$x.orig[int.ind],
               s.x = x + curve[[pa]]$s.x[int.ind],
               #lower = x + curve[[pa]]$lower[int.ind],
               #upper = x + curve[[pa]]$upper[int.ind],
               type = ind);})
    # Label as treatment group
    ctrl <- FALSE
  }

  # Combine into one dataframe for each pa
  curve.with.int[[pa]] <- do.call("rbind", curve.with.int[[pa]])
}

# Combine into one list of dataframes for each PA variable
final.curve <- list()
final.curve[[pa_vars[1]]] <- do.call("rbind", lapply(names(curve.with.int)[which(grepl("MVPA", names(curve.with.int)))],
                                                     function(x) curve.with.int[[x]]))
final.curve[[pa_vars[2]]] <- do.call("rbind", lapply(names(curve.with.int)[which(grepl("ABD", names(curve.with.int)))],
                                                     function(x) curve.with.int[[x]]))

####################################################################
# Plot the Additive Terms
####################################################################

## Plot the estimated curves
glist <- list()
J <- 1 # loop count
for(pa in pa_vars){

  # Get x-axis label for plotting
  xlab <- paste0("Difference in ", plot.title[J], " (min/day)")

  # Get y-axis label for plotting
  if(grepl("intensity", y)){
    ylab <- "Difference in fatigue intensity after baseline"
  }
  else{
    ylab <- "Difference in fatigue interference after baseline"
  }
  # Organize the data frame for plotting by order of x-axis
  est.df <- final.curve[[pa]][order(final.curve[[pa]][,1]),]
  est.df <- distinct(est.df)

  # Create levels of each curve to specify group type and time period
  est.df$type <- as.factor(est.df$type)
  levels(est.df$type) <- c("Ctrl-M3", "Ctrl-M6", "Ctrl-M12", "Trt-M3", "Trt-M6", "Trt-M12")

  # Create plot and save to list
  glist[[pa]] <- ggplot(est.df, aes(x = x, y = (s.x), group = type, color = type)) +
    geom_line() + theme_bw() +
    geom_rug(aes(x = x), sides = "b", color = "black") +
    labs(x = xlab, y = ylab, title = plot.title[J]) +
    theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
          text = element_text(size = 10), plot.subtitle = element_text(hjust = 0.5))

  # Update loop count
  J <- J + 1
}

# pdf(file = "./figures/factor_by_curve_varying_intercept/EstimatedCurve_Intensity_JAGS.pdf",
#     width = 10, height = 4.5)
png(file="./plots/2var/fatigue_comparison.png", width=600, height=350)
do.call("grid.arrange", c(glist, nrow = floor(sqrt(length(glist)))))
dev.off()

####################################################################
## MCMC diagnostics
####################################################################

## Obtain posterior samples as a mcmc.list object
samples.coda <- coda.samples(jags.mod, c("beta", "g1 ", "g2","g3","g4","sigma_u1","sigma_u2",
                                         "sigma_u3","sigma_u4","sigma_y"), n.iter = 50000, thin = 50)
samples.array <- as.array(samples.coda)
dimnames(samples.array)[[2]] <- c(colnames(X.mat),"MVPA_C", "MVPA_I","ABD_C","ABD_I","sigma_u1 for MVPA_C", "sigma_u2 for MVPA_I","sigma_u3 for ABD_C", "sigma_u4 for ABD_I", "sigma_y")

##matrix to store the g estimates
g_mat = matrix(NA, 1000,4)

for(i in 1:4){
  g_mat[,i ] <- samples.array[ ,(ncol(X.mat) +i)]

}

##creating the trace plots of the g estimates and saving them
gamma.list <- lapply(1:4, function(i) {
  ggplot(data = as.data.frame(g_mat[ ,i]), aes(x = 1:1000, y = g_mat[ ,i] ))+
    geom_line() + xlab("Index") +  theme_bw()+
    ggtitle(colnames(samples.array)[(ncol(X.mat) +i)])+ylab(NULL)
})
png(file="./plots/2var/gamma.png")
ggarrange(plotlist = gamma.list, nrow = 2, ncol = 2)
dev.off()

##matrix to store the beta estimates
beta_mat = matrix(NA, 1000, ncol(X.mat))
for(i in 1:ncol(X.mat)){
  beta_mat[ , i] =  samples.array[ ,i]
}
##creating the trace plots of the beta estimates and saving them
beta.list <- lapply(1:ncol(X.mat), function(i) {
  ggplot(data = as.data.frame(beta_mat[ ,i]), aes(x = 1:1000, y = beta_mat[ ,i] ))+
    geom_line() + xlab("Index") +  theme_bw()+
    ggtitle(colnames(samples.array)[i])+ylab(NULL)
})

png(file="./plots/2var/beta.png")
ggarrange(plotlist = beta.list, nrow = 4, ncol = 5)
dev.off()

##creating the trace plots of the sigma estimates and saving them
sigmas.list <- lapply((ncol(X.mat)+5):ncol(samples.array), function(i) {
  ggplot(data = as.data.frame(samples.array[ ,i]), aes(x = 1:1000, y = samples.array[ ,i] ))+
    geom_line() + xlab("Index")+  theme_bw() +
    ggtitle(colnames(samples.array)[i])+ylab(NULL)
})
png(file="./plots/2var/sigma_plots.png")
ggarrange(plotlist = sigmas.list, nrow = 3, ncol = 2)
dev.off()



##########################################################################
# END
##########################################################################





## Organize the linear parameter results
# control = c(1,3:7,8,10,11,14:16)
# trt = c(2,3:7,9,12,13,17:19)
#
#
#
# est_1 = rep(apply(X.mat,2,mean)%*% apply(samples.jags$beta, 1, mean),length(int1)) + MVPA_C[int1]*rep(apply(samples.jags$g1, 1, mean),length(int1)) +
#   rep(mean(MVPA_I)*apply(samples.jags$g2, 1, mean),length(int1)) + rep(mean(ABD_C)*apply(samples.jags$g3, 1, mean),length(int1)) + rep(mean(ABD_I)*apply(samples.jags$g4, 1, mean),length(int1))
# est_2 = rep(apply(X.mat,2,mean)%*% apply(samples.jags$beta, 1, mean),length(int2)) + MVPA_C[int2]*rep(apply(samples.jags$g1, 1, mean),length(int2)) +
#   rep(mean(MVPA_I)*apply(samples.jags$g2, 1, mean),length(int2)) + rep(mean(ABD_C)*apply(samples.jags$g3, 1, mean),length(int2)) + rep(mean(ABD_I)*apply(samples.jags$g4, 1, mean),length(int2))
# est_3 = rep(apply(X.mat,2,mean)%*% apply(samples.jags$beta, 1, mean),length(int3)) + MVPA_C[int3]*rep(apply(samples.jags$g1, 1, mean),length(int3)) +
#   rep(mean(MVPA_I)*apply(samples.jags$g2, 1, mean),length(int3)) + rep(mean(ABD_C)*apply(samples.jags$g3, 1, mean),length(int3)) + rep(mean(ABD_I)*apply(samples.jags$g4, 1, mean),length(int3))
#
# x = MVPA_C[c(int1, int2, int3)]
# y.est = c(est_1,est_2,est_3)
# type = as.factor(rep(c("Ctrl-M3", "Ctrl-M6", "Ctrl-M12"), c(length(int1), length(int2),length(int3))))
# plot.data = data.frame(x,y.est,type)
#
#
# ggplot(plot.data, aes(x = x)) + geom_line(aes(y = est_1)) + geom_line(aes(y = est_2))
#
# ggplot(plot.data, aes(x = 1:148, y = y.est, group = type, color = type)) +
#   geom_line() + theme_bw()
# +
#   geom_rug(aes(x = x), sides = "b", color = "black") +
#   labs(x = xlab, y = ylab, title = plot.title[J]) +
#   theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
#         text = element_text(size = 10), plot.subtitle = element_text(hjust = 0.5))


