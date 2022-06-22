dataPreprocess = function(sim_data, type=NULL, knots=NULL, basis=NULL, nci=NULL,control.points=NULL ){

  if(type == "2var"){
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

    return(list("X.mat" = X.mat, "MVPA_C" = MVPA_C, "MVPA_I" = MVPA_I, "ABD_C" = ABD_C, "ABD_I" = ABD_I))
  }else{
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


  # Formula for model for mean ABD and MVPA on fatigue
  k <- knots # number of knots for spline
  form <- paste0(y_sc, " ~ -1 + int.time")
  for(pa in pa_vars){
    form <- paste0(form, " + s(", pa, ", bs = basis, k = ", k, ", by = as.factor(Randomization))")
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
  X.mat <- X.mat.full[ ,1:nci] #covariates
  W1_C.mat <- X.mat.full[ ,(nci+1):(nci+control.points-1)]  # MVPA for control group for different knot points
  W1_I.mat <- X.mat.full[ ,(nci+control.points):(nci+2*control.points-2)]  # MVPA for treatment group for different knot points
  W2_C.mat <- X.mat.full[ ,(nci+2*control.points-1):(nci+3*control.points-3)]  # ABD for control group for different knot points
  W2_I.mat <- X.mat.full[ ,(nci+3*control.points-2):(nci+4*control.points-4)]  # ABD for treatment group for different knot points

  S1 <- cbind(fit$smooth[[1]]$S[[1]], fit$smooth[[1]]$S[[2]]) ## penalty matrices for MVPA control group
  S2 <- cbind(fit$smooth[[2]]$S[[1]], fit$smooth[[2]]$S[[2]]) ## penalty matrices for MVPA treatment group
  S3 <- cbind(fit$smooth[[3]]$S[[1]], fit$smooth[[3]]$S[[2]]) ## penalty matrices for ABD control group
  S4 <- cbind(fit$smooth[[4]]$S[[1]], fit$smooth[[4]]$S[[2]]) ## penalty matrices for ABD treatment group

  return(list("X.mat"=X.mat,"W1_C.mat"=W1_C.mat,"W1_I.mat"=W1_I.mat, "W2_C.mat"=W2_C.mat,"W2_I.mat"=W2_I.mat,
              "S1"=S1,"S2"=S2,"S3"=S3,"S4"=S4))
}
}
