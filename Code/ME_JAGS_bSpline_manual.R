################################################################################
# This program creates a measurement error model on the simulated
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
library(ggpubr)
library(splines)
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

    temp1[i] ~ dnorm(0, taue1)
    temp2[i] ~ dnorm(0, taue2)
    temp3[i] ~ dnorm(0, taue3)
    temp4[i] ~ dnorm(0, taue4)
  }
  # S(MVPA) model for control group
  Z1_C = bs(temp1, knots = (k-2))

  # S(MVPA) model for treatment group
  Z1_I = bs(temp2, knots = (k-2))

  # S(mean ABD) model for  control group
  Z2_C = bs(temp3, knots = (k-2))

  # S(mean ABD) model for treatment group
  Z2_I = bs(temp4, knots = (k-2))


 # classical measurement error model
 for(i in 1:n){
  for(j in 1:(k-1)){
    W1_C[i,j] ~ dnorm(z1_C[i, j], tauu1)
    W1_I[i,j] ~ dnorm(z1_I[i, j], tauu2)
    W2_C[i,j] ~ dnorm(z2_C[i, j] , tauu3)
    W2_I[i,j] ~ dnorm(z2_I[i, j], tauu4)

  }
 }

  ## Prior for covariates and intercept
  for (i in 1:(nc+ni)) { b[i] ~ dnorm(0,0.0075) }

  ## Prior for spline parameters
  for(i in 1:(k-1)){
  g1[i] ~ dnorm(0,0.01)
  g2[i] ~ dnorm(0,0.01)
  g3[i] ~ dnorm(0,0.01)
  g4[i] ~ dnorm(0,0.01)
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

  sigma_u1 <- 1/tauu1
  sigma_u2 <- 1/tauu2
  sigma_u3 <- 1/tauu3
  sigma_u4 <- 1/tauu4
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
                  nc = nc, k = k, X = X.mat, W1_C = W1_C.mat, W1_I = W1_I.mat, W2_C = W2_C.mat,
                  W2_I = W2_I.mat, zero = rep(0, ncol(X.mat.full)))

# Intilialize JAGS model
jm <- jags.model(textConnection(model_string), data = data.jags,
                 inits = list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 2021),
                 n.chains = 1, n.adapt = 100)

# specify number of burn-in
update(jm, n.burn = 200)

# Obtain JAGS samples
start = Sys.time()
samples <- jags.samples(jm, c("b", "g1", "g2", "g3", "g4","sigma_u1","sigma_u2",
                              "sigma_u3","sigma_u4","scale"), n.iter = 1000, thin = 2)
stop = Sys.time()
time.elaspsed = stop - start; time.elaspsed

samples.jags <- c(apply(samples$b, 1, mean),apply(samples$g1, 1, mean),apply(samples$g2, 1, mean),
                  apply(samples$g3, 1, mean),apply(samples$g4, 1, mean))
####################################################################
## MCMC diagnostics
####################################################################

# ## Obtain posterior samples as a mcmc.list object
# samples.coda <- coda.samples(jm, c("b", "g1", "g2", "g3", "g4","sigma_u1","sigma_u2",
#                                    "sigma_u3","sigma_u4","scale"), n.iter = 50000, thin = 50)
samples.array <- data.frame(t(samples$b[ , ,1]), t(samples$g1[ , ,1]), t(samples$g2[ , ,1]), t(samples$g3[ , ,1]), t(samples$g4[ , ,1]),
                            samples$sigma_u1[ , ,1], samples$sigma_u2[ , ,1], samples$sigma_u3[ , ,1],
                            samples$sigma_u4[ , ,1], samples$scale[ , ,1])
for(i in 1:ncol(X.mat)){
  colnames(samples.array)[i] <- colnames(X.mat)[i]
}
for(i in 1:10){
  colnames(samples.array)[(ncol(X.mat)+i)] <- paste("MVPA_C_k_",i)
  colnames(samples.array)[(ncol(X.mat)+(k-1)+i)] <- paste("MVPA_I_k_",i)
  colnames(samples.array)[(ncol(X.mat)+2*(k-1)+i)] <- paste("ABD_C_k_",i)
  colnames(samples.array)[(ncol(X.mat)+3*(k-1)+i)] <- paste("ABD_I_k_",i)
}
colnames(samples.array)[61:65] <- c("sigma_u1","sigma_u2", "sigma_u3","sigma_u4","sigma_y")

write.csv(samples.array, "C:/Users/anany/Desktop/JAGSME/Code/samples.bspline.csv")


#matrices to store the g values for the (k-1) knots
n_row = ncol(samples$g1)
g1_mat = matrix(NA, n_row,10)
g2_mat = matrix(NA, n_row,10)
g3_mat = matrix(NA, n_row,10)
g4_mat = matrix(NA, n_row,10)
for(i in 1:10){
  g1_mat[,i ] <- samples$g1[i, ,1]
  g2_mat[,i ] <- samples$g2[i, ,1]
  g3_mat[,i ] <- samples$g3[i, ,1]
  g4_mat[,i ] <- samples$g4[i, ,1]
}

#creating the list of traceplots corresponding to the posterior distributions
#at each knot point for each g parameter
g1_list <- lapply(1:10, function(i) {
  ggplot(data = as.data.frame(g1_mat[ ,i]), aes(x = 1:n_row, y = g1_mat[ ,i] ))+
    geom_line() + xlab("Index (B-S)")+
    ggtitle(paste("MVPA_C_k_",i))+ylab(NULL)
})

g2_list <- lapply(1:10, function(i) {
  ggplot(data = as.data.frame(g2_mat[ ,i]), aes(x = 1:n_row, y = g2_mat[ ,i] ))+
    geom_line() + xlab("Index (B-S)")+
    ggtitle(paste("MVPA_I_k_",i))+ylab(NULL)
})
g3_list <- lapply(1:10, function(i) {
  ggplot(data = as.data.frame(g3_mat[ ,i]), aes(x = 1:n_row, y = g3_mat[ ,i] ))+
    geom_line() + xlab("Index (B-S)")+
    ggtitle(paste("ABD_C_k_",i))+ylab(NULL)
})
g4_list <- lapply(1:10, function(i) {
  ggplot(data = as.data.frame(g4_mat[ ,i]), aes(x = 1:n_row, y = g4_mat[ ,i] ))+
    geom_line() + xlab("Index (B-S)")+
    ggtitle(paste("ABD_I_k_",i))+ylab(NULL)
})

#saving the plots as .png files
#png(file="./plots/Bsplines/MVPA_C.png")
ggarrange(plotlist = g1_list, nrow = 3, ncol = 4)
#dev.off()
#png(file="./plots/Bsplines/MVPA_I.png")
ggarrange(plotlist = g2_list, nrow = 3, ncol = 4)
#dev.off()
#png(file="./plots/Bsplines/ABD_C.png")
ggarrange(plotlist = g3_list, nrow = 3, ncol = 4)
#dev.off()
#png(file="./plots/Bsplines/ABD_I.png")
ggarrange(plotlist = g4_list, nrow = 3, ncol = 4)
#dev.off()

#creating the list of traceplots corresponding to the posterior distributions
#of each beta parameter ordered by interval time, age and BMI, hormone treatment
# and cancer sample

beta_mat = matrix(NA, n_row, 20)
for(i in 1:20){
  beta_mat[ , i] =  samples$b[i, ,1]
}

colnames(beta_mat) = colnames(X.mat)
covariates <- colnames(X.mat)

int.time.list <- lapply(1:6, function(i) {
  ggplot(data = as.data.frame(beta_mat[ ,i]), aes(x = 1:n_row, y = beta_mat[ ,i] ))+
    geom_line() + xlab("Index (B-S)")+
    ggtitle(covariates[i])+ylab(NULL)
})

ageBMI.list <- lapply(c(7,8,14,15), function(i) {
  ggplot(data = as.data.frame(beta_mat[ ,i]), aes(x = 1:n_row, y = beta_mat[ ,i] ))+
    geom_line() + xlab("Index (B-S)")+
    ggtitle(covariates[i])+ylab(NULL)
})

ht_sample.list <- lapply(c(9,10,16,17), function(i) {
  ggplot(data = as.data.frame(beta_mat[ ,i]), aes(x = 1:n_row, y = beta_mat[ ,i] ))+
    geom_line() + xlab("Index (B-S)")+
    ggtitle(covariates[i])+ylab(NULL)
})

breast_sample.list <- lapply(c(11:13,18:20), function(i) {
  ggplot(data = as.data.frame(beta_mat[ ,i]), aes(x = 1:n_row, y = beta_mat[ ,i] ))+
    geom_line() + xlab("Index (B-S)")+
    ggtitle(covariates[i])+ylab(NULL)
})

#saving the plots as .png files
#png(file="./plots/Bsplines/int.time.png")
ggarrange(plotlist = int.time.list, nrow = 3, ncol = 2)
#dev.off()
#png(file="./plots/Bsplines/ageBMI.png")
ggarrange(plotlist = ageBMI.list, nrow = 2, ncol = 2)
#dev.off()
#png(file="./plots/Bsplines/ht_sample.png")
ggarrange(plotlist = ht_sample.list, nrow = 2, ncol = 2)
#dev.off()
#png(file="./plots/Bsplines/breast_sample.png")
ggarrange(plotlist = breast_sample.list, nrow = 3, ncol = 2)
##dev.off()

#creating the list of traceplots corresponding to the posterior distributions
#of the precision parameters of the overall fatigue intensity and the
#measurement error variable
sigma_mat = matrix(NA, n_row, 5)
sigma_mat[ ,1] <- samples$sigma_u1[1, ,1]
sigma_mat[ ,2] <- samples$sigma_u2[1, ,1]
sigma_mat[ ,3] <- samples$sigma_u3[1, ,1]
sigma_mat[ ,4] <- samples$sigma_u4[1, ,1]
sigma_mat[ ,5] <- samples$scale[1, ,1]
s2_names <- c("sigma_u1 for MVPA_C", "sigma_u2 for MVPA_I","sigma_u3 for ABD_C", "sigma_u4 for ABD_I", "sigma_y")

sigmas.list <- lapply(1:5, function(i) {
  ggplot(data = as.data.frame(sigma_mat[ ,i]), aes(x = 1:n_row, y = sigma_mat[ ,i] ))+
    geom_line() + xlab("Index(B-S)")+  theme_bw() +
    ggtitle(s2_names[i])+ylab(NULL)
})

#saving the plots as .png files
#png(file="./plots/Bsplines/sigma_plots.png")
ggarrange(plotlist = sigmas.list, nrow = 2, ncol = 3)
##dev.off()
