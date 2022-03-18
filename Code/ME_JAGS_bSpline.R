################################################################################
# This progrma creates a factor-by-curve interaction model on the simulated
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
library(bSpline)
# Source functions
function_directory <- getwd()
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
zC <- as.matrix(dataC[ , covars[grepl(".C", covars)]])
zI <- as.matrix(dataI[ , covars[grepl(".I", covars)]])
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
# Define JAGS model
#################################################################################
## nc: number of covariates
## ni: number of intercepts
## k: number of knots per spline

model_string <- "model {

  # Likelihood (deign matrix %*% matrix of parameters to estimate)
  mu <- X %*% b + z1_C %*% g1 + z1_I %*% g2 + z2_C %*% g3 + z2_I %*% g4  ## expected response
                                                                         ## Z1 = actual MVPA for trt & control groups
                                                                         ## Z2 = actual ABD for trt & control groups
          ### temp's are the values of the basis functions evaluated at that Z

  for (i in 1:n){
  ## response model
  y[i] ~ dnorm(mu[i],tau)


  for(j in 1:(k-1)){
  # S(MVPA) model for control group
  z1_C[i, j] ~ dnorm(0, taue1)
  # classical measurement error model
  W1_C[i,j] ~ dnorm(z1_C[i, j], tauu1)

  # S(MVPA) model for treatment group
  z1_I[i, j] ~ dnorm(0, taue2)
  # classical measurement error model
  W1_I[i,j] ~ dnorm(z1_I[i, j], tauu2)


  # s(mean ABD) model for  control group
  z2_C[i, j] ~ dnorm(0, taue3)
  # classical measurement error model
  W2_C[i,j] ~ dnorm(z2_C[i, j] , tauu3)


  # s(mean ABD) model for treatment group
  z2_I[i, j] ~ dnorm(0, taue4)
  # classical measurement error model

  W2_I[i,j] ~ dnorm(z2_I[i, j], tauu4)

  }

  }


  ## Prior for covariates and intercept
  for (i in 1:(nc+ni)) { b[i] ~ dnorm(0,0.0075) }

  ## Prior for s(MVPA) control group
  K1 <- S1[1:(k-1),1:(k-1)] * lambda[1]  + S1[1:(k-1),k:(2*(k-1))] * lambda[2]
  g1[1:(k-1)] ~ dmnorm(zero[(nc+ni+1):(nc+ni+k-1)],K1)

  ## Prior for s(MVPA) treatment group
  K2 <- S2[1:(k-1),1:(k-1)] * lambda[1]  + S2[1:(k-1),k:(2*(k-1))] * lambda[2]
  g2[1:(k-1)] ~ dmnorm(zero[(nc+ni+k):(nc+ni+2*(k-1))],K2)

  ## Prior for s(mean ABD) control group
  K3 <- S3[1:(k-1),1:(k-1)] * lambda[1]  + S3[1:(k-1),k:(2*(k-1))] * lambda[2]
  g3[1:(k-1)] ~ dmnorm(zero[(nc+ni+2*(k-1)+1):(nc+ni+3*(k-1))],K3)

  ## Prior for s(mean ABD) treatment group
  K4 <- S4[1:(k-1),1:(k-1)] * lambda[1]  + S4[1:(k-1),k:(2*(k-1))] * lambda[2]
  g4[1:(k-1)] ~ dmnorm(zero[(nc+ni+3*(k-1)+1):(nc+ni+4*(k-1))],K4)

  ## Smoothing parameter priors
  for (i in 1:2) {
    lambda[i] ~ dgamma(1, 1)
    rho[i] <- log(lambda[i])
  }

  tau ~ dgamma(.05,.005) ## precision parameter prior
  scale <- 1/tau ## convert tau to standard GLM scale

  # prior distributions error model
  tauu1 ~ dgamma(0.01, 0.01)
  tauu2 ~ dgamma(0.01, 0.01)
  tauu3 ~ dgamma(0.01, 0.01)
  tauu4 ~ dgamma(0.01, 0.01)
  taue1 ~  dgamma(0.01, 0.01)
  taue2 ~  dgamma(0.01, 0.01)
  taue3 ~  dgamma(0.01, 0.01)
  taue4 ~  dgamma(0.01, 0.01)
}"

#################################################################################
# Build the JAGS model
#################################################################################

# Formula for model for mean ABD and MVPA on fatigue
k <- 11 # number of knots for spline
form <- paste0(y_sc, " ~ -1 + int.time")
for(pa in pa_vars){
  form <- paste0(form, " + s(", pa, ", bs = 'bs', k = ", k, ", by = as.factor(Randomization))")
}

for(param in c(colnames(zC), colnames(zI))){
  form <- paste0(form, " + ", param)
}

# Fit model for mean ABD and MVPA on fatigue intensity
fit <- gam(as.formula(form), data = data.sim,
           method = 'REML', select = TRUE)

jags.file <- paste(function_directory, "/code/","/test.jags",sep="")
fit.jagam <- jagam(as.formula(form),data=data.sim, file=jags.file,
                   sp.prior="gamma",diagonalize=TRUE)
# Derive linear predictor matrix
X.mat.full <- predict.gam(fit, type = "lpmatrix")
X.mat <- X.mat.full[ ,1:20] #covariates
W1_C.mat <- X.mat.full[ ,21:30]  # MVPA for control group for different knot points
W1_I.mat <- X.mat.full[ ,31:40]  # MVPA for treatment group for different knot points
W2_C.mat <- X.mat.full[ ,41:50]  # ABD for control group for different knot points
W2_I.mat <- X.mat.full[ ,51:60]  # ABD for treatment group for different knot points

S1 <- cbind(fit$smooth[[1]]$S[[1]], fit$smooth[[1]]$S[[2]]) ## penalty matrices for MVPA control group
S2 <- cbind(fit$smooth[[2]]$S[[1]], fit$smooth[[2]]$S[[2]]) ## penalty matrices for MVPA treatment group
S3 <- cbind(fit$smooth[[3]]$S[[1]], fit$smooth[[3]]$S[[2]]) ## penalty matrices for ABD control group
S4 <- cbind(fit$smooth[[4]]$S[[1]], fit$smooth[[4]]$S[[2]]) ## penalty matrices for ABD treatment group


## Build the final dataset to fit Bayesian model
ni <- length(unique(data$int.time)) # number of intercept terms
nc <- ncol(X.mat.full) - 4*(k-1) - ni # number of covariates
data.jags <- list(y = as.vector(Y), n = length(Y), ni = ni,
                  nc = nc, k = k, X = X.mat, W1_C = W1_C.mat, W1_I = W1_I.mat, W2_C = W2_C.mat, W2_I = W2_I.mat, S1 = S1,
                  S2 = S2, S3 = S3, S4 = S4, zero = rep(0, ncol(X.mat.full)))

# Intilialize JAGS model
jm <- jags.model(textConnection(model_string), data = data.jags,
                 inits = list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 2021),
                 n.chains = 1, n.adapt = 1000)

# specify number of burn-in
update(jm, n.burn = 5000)

# Obtain JAGS samples
samples <- jags.samples(jm, c("b", "g1", "g2", "g3", "g4","scale"), n.iter = 3000, thin = 3)



####################################################################
## MCMC diagnostics
####################################################################

## Obtain posterior samples as a mcmc.list object
samples.coda <- coda.samples(jm, c("b", "g1", "g2", "g3", "g4","scale"), n.iter = 50000, thin = 50)
samples.array <- as.array(samples.coda)
dimnames(samples.array)[[2]][1:20] <- colnames(X.mat)

# ## Traceplots for covariates
# par(mfrow = c(4, 5))
# traceplot(mcmc(samples.array[,ind.cov]))

## Traceplots for spline parameters
par(mfrow = c(3,3))
for(j in 1:61){
  traceplot(mcmc(samples.array[,j]),main = dimnames(samples.array)$var[j])
}


## effect sample sizes
par(mfrow=c(1,1))
Neff <- effectiveSize(samples.coda)
plot(1:length(Neff), Neff, xlab = "parameter")

################################################################################
# Or create JAGS code automatically to fit model using jagam()
################################################################################
# # Create JAGS code for model
# jd <- jagam(as.formula(form), data = data.sim, family = gaussian,
#             file = "./code/FBC_VI_model_test.jags")
# jd$jags.data$y <- as.vector(jd$jags.data$y)
# data.jags <- jd$jags.data
#
# # Initialize JAGS model
# jm <- jags.model("./code/FBC_VI_model_test.jags", data = data.jags,
#                  inits = list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 2021),
#                  n.chains = 1, n.adapt = 1000)
#
# # specify number of burn-in
# update(jm, n.burn = 5000)
#
# # Obtain samples
# samples <- jags.samples(jm, c("b"), n.iter = 3000, thin = 3)

################################################################################
# Organize the results
################################################################################

## Organize the linear parameter results

# Obtain index of covariates and intercept
ind.cov <- 1:(ni+nc)

# Get mean and sd from MCMC samples
theta.sample <- samples$b[ind.cov,,1]
est <- apply(theta.sample, 1, mean)
se <- apply(theta.sample, 1, sd)

# Calculate p-value for whether estimate != 0
pvals <-  2 * (1 - pnorm(abs(est / se)))

## Combine output and print results
results <- cbind(est, se, pvals)
rownames(results) <- colnames(X.mat)[ind.cov]
print(results, digits = 2)


## Organize the curve output

# Combine scaled PA vars into matrix
X <- data.sim[ , pa_vars]

# Combine original units of PA vars into matrix
X.orig <- rbind(dataC[c(t1.C, t2.C, t3.C), pa_vars_orig],
                dataI[c(t1.I, t2.I, t3.I), pa_vars_orig])
# Number of spline bases per term
nb <- k - 1

# List to save curve results
curve <- list()

# Loop over each smooth function
# (one curve for each group and for each PA variable)
for(i in 1:(2*length(pa_vars))){

  # Data is grouped as MVPA control group, MVPA trt group, ABD control group, ABD trt group
  # Specify type = 0 for control in data groups 1, 3
  if(i %in% c(1, 3)){
    type = 0
  }
  # Specify type = 1 for treatment in data groups 2, 4
  else{
    type = 1
  }
  # Indices of corresponding spline basis in X.spl
  ind.spl <- nc + ni + ((i-1)*nb+1):(i*nb)

  # Calculate estimated curve for each MCMC sample
  fHat.all <- X.mat.full[,ind.spl] %*% apply(samples$g1, 1, mean)

  # Backtransform estimated curves into original units
  fHat.tilde <- fHat.all * sd(data[ , y]) + mean(data[ , y])

  # Obtain mean and 95% CI for curve
  fHat.mean <- apply(fHat.tilde, 1, mean)
  lower <- apply(fHat.tilde, 1, function(x) quantile(x, 0.025))
  upper <- apply(fHat.tilde, 1, function(x) quantile(x, 0.975))


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
  curve[[paste0(pa, i)]] <- data.frame(x.orig = x.orig, s.x = fHat.mean)

  # Smooth CI estimates and save into list
  curve[[paste0(pa, i)]]$lower <- gam(lower ~ s(x), method = "REML")$fitted.values
  curve[[paste0(pa, i)]]$upper <- gam(upper ~ s(x), method = "REML")$fitted.values

}

# Obtain intercepts in original scale
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
               lower = x + curve[[pa]]$lower[int.ind],
               upper = x + curve[[pa]]$upper[int.ind],
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
               lower = x + curve[[pa]]$lower[int.ind],
               upper = x + curve[[pa]]$upper[int.ind],
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
  glist[[pa]] <- ggplot(est.df, aes(x = x, y = s.x, group = type, color = type)) +
    geom_line() + theme_bw() + ylim(-1.8, 2.6) +
    geom_rug(aes(x = x), sides = "b", color = "black") +
    labs(x = xlab, y = ylab, title = plot.title[J]) +
    theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
          text = element_text(size = 10), plot.subtitle = element_text(hjust = 0.5))

  # Update loop count
  J <- J + 1
}

# pdf(file = "./figures/factor_by_curve_varying_intercept/EstimatedCurve_Intensity_JAGS.pdf",
#     width = 10, height = 4.5)
png(file="./plots/Bsplines/fatigue_comparison.png", width=600, height=350)
do.call("grid.arrange", c(glist, nrow = floor(sqrt(length(glist)))))
dev.off()


g1_mat = matrix(NA, 1000,10)
g2_mat = matrix(NA, 1000,10)
g3_mat = matrix(NA, 1000,10)
g4_mat = matrix(NA, 1000,10)
for(i in 1:10){
  g1_mat[,i ] <- samples.array[ ,20+i]
  g2_mat[,i ] <- samples.array[ ,30+i]
  g3_mat[,i ] <- samples.array[ ,40+i]
  g4_mat[,i ] <- samples.array[ ,50+i]
}
colnames(g1_mat) <- colnames(W1_C.mat)
colnames(g2_mat) <- colnames(W1_I.mat)
colnames(g3_mat) <- colnames(W2_C.mat)
colnames(g4_mat) <- colnames(W2_I.mat)

g1_list <- list()
g2_list <- list()
g3_list <- list()
g4_list <- list()
for(cov in 1:10){
  g1_list[[cov]] <- ggplot(as.data.frame(g1_mat[ ,cov]), aes(x = 1:1000, y = g1_mat[ ,cov])) + geom_line() + xlab("Index")+
    ggtitle(paste("MVPA_C_knot",cov))+ylab("")
  g2_list[[cov]] <- ggplot(as.data.frame(g2_mat[ ,cov]), aes(x = 1:1000, y = g2_mat[ ,cov])) + geom_line() + xlab("Index")+
    ggtitle(paste("MVPA_I_knot",cov))+ylab("")
  g3_list[[cov]] <- ggplot(as.data.frame(g3_mat[ ,cov]), aes(x = 1:1000, y = g3_mat[ ,cov])) + geom_line() + xlab("Index")+
    ggtitle(paste("ABD_C_knot",cov))+ylab("")
  g4_list[[cov]] <- ggplot(as.data.frame(g4_mat[ ,cov]), aes(x = 1:1000, y = g4_mat[ ,cov])) + geom_line() + xlab("Index")+
    ggtitle(paste("ABD_I_knot",cov))+ylab("")
}
png(file="./plots/Bsplines/MVPA_C.png", width=600, height=350)
do.call("grid.arrange", c(g1_list, nrow = 3))
dev.off()
png(file="./plots/Bsplines/MVPA_I.png", width=600, height=350)
do.call("grid.arrange", c(g2_list, nrow = 3))
dev.off()
png(file="./plots/Bsplines/ABD_C.png", width=600, height=350)
do.call("grid.arrange", c(g3_list, nrow = 3))
dev.off()
png(file="./plots/Bsplines/ABD_I.png", width=600, height=350)
do.call("grid.arrange", c(g4_list, nrow = 3))
dev.off()

beta_mat = matrix(NA, 1000, 20)
for(i in 1:20){
  beta_mat[ , i] =  samples.array[ ,i]
}

colnames(beta_mat) = colnames(X.mat)
int.time.list <- list()
ageBMI.list <- list()
ht_sample.list <- list()
breast_sample.list <- list()


covariates <- colnames(X.mat)
for(cov in 1:6){
  int.time.list[[cov]] <- ggplot(as.data.frame(beta_mat[ ,cov]), aes(x = 1:1000, y = beta_mat[ ,cov])) + geom_line() + xlab("Index")+
    ggtitle(paste(covariates[cov]))+ylab("")
}

k = 1
for(cov in c(7,8,14,15)){
  #print(paste("cov = ", cov, head(beta_mat[ ,cov])))
  ageBMI.list[[k]] <- ggplot(as.data.frame(beta_mat[ ,cov]), aes(x = 1:1000, y = beta_mat[ ,cov])) + geom_line() + xlab("Index")+
    ggtitle(paste(covariates[cov]))+ylab("")
  k = k+1
}

k = 1
for(cov in c(9,10,16,17)){
  ht_sample.list[[k]] <- ggplot(as.data.frame(beta_mat[ ,cov]), aes(x = 1:1000, y = beta_mat[ ,cov])) + geom_line() + xlab("Index")+
    ggtitle(paste(covariates[cov]))+ylab("")
  k = k+1
}

k = 1
for(cov in c(11:13, 18:20)){
  breast_sample.list[[k]] <- ggplot(as.data.frame(beta_mat[ ,cov]), aes(x = 1:1000, y = beta_mat[ ,cov])) + geom_line() + xlab("Index")+
    ggtitle(paste(covariates[cov]))+ylab("")
  k = k+1
}
png(file="./plots/Bsplines/int.time.png", width=600, height=350)
do.call("grid.arrange", c(int.time.list, nrow = 3))
dev.off()
png(file="./plots/Bsplines/ageBMIfatigue_comparison.png", width=600, height=350)
do.call("grid.arrange", c(ageBMI.list, nrow = 2))
dev.off()
png(file="./plots/Bsplines/ht_sample.png", width=600, height=350)
do.call("grid.arrange", c(ht_sample.list, nrow = 2))
dev.off()
png(file="./plots/Bsplines/breast_sample.png", width=600, height=350)
do.call("grid.arrange", c(breast_sample.list, nrow = 3))
dev.off()
