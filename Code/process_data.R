################################################################################
# This program pre-processes the data into analysis-ready format
################################################################################

# data: dataframe to analyze
# vars: vector of variable names (response and PA variables) to take difference from baseline
# group: variable indicating group assignment

preprocess_data <- function(data, vars, group){
  
  ## Create variables for difference in response from baseline
  
  # Subset data without baseline observations
  data.after.base <- data[which(data$svy_interval != 0), ]
  
  # Loop over each subject
  for(i in unique(data$new_id)){
    # Subset data by subject
    data.subj <- data[data$new_id == i, ]
    
    # If there exists a baseline measure
    if(any(data.subj$svy_interval == 0)){
      # Take difference from baseline
      diff <- apply(data.subj[which(data.subj$svy_interval != 0), vars], 1, 
                    function(x) x - data.subj[which(data.subj$svy_interval == 0), vars])
    }
    # If there is not a baseline measure
    else{
      diff <- c()
    }
    
    
    # If there is more than just a baseline observation, add to vector
    if(length(diff) != 0){
      # Combine into matrix 
      diff.mat <- do.call(rbind, diff)
      colnames(diff.mat) <- paste0(vars, "_diff")
      
      # Add to data frame
      data.after.base[which(data.after.base$new_id == i), colnames(diff.mat)] <- diff.mat
    }
  }
  
  # The subjects that are missing baseline measurements
  id.no.base <- unique(data.after.base$new_id[which(is.na(data.after.base[ , 'fatigue_mean_intensity_score_diff']))])
  
  # Discard any NA values in new data
  data.after.base <- data.after.base[-which(data.after.base$new_id %in% id.no.base), ]
  
  ## Scale the response and PA measures
  for(variable in colnames(diff.mat)){
    data.after.base[ , paste0(variable, "_sc")] <- scale(data.after.base[ , variable]) 
    # var.mat <- cbind(scale(ifelse(data.after.base[, group] == 0, data.after.base[ , variable], NA)),
    #                  scale(ifelse(data.after.base[, group] == 1, data.after.base[ , variable], NA)))
    # data.after.base[ , paste0(variable, "_sc")] <- rowSums(var.mat, na.rm = T)
  }
  
  # Create intercept indicator for time period and group
  # 1-3: control group for time periods M3, M6, M12
  # 4-5: treatment group for time periods M3, M6, M12
  data.after.base$int.time <- rep(1, length = nrow(data.after.base))
  for(i in 1:nrow(data.after.base)){
    if(data.after.base[, group][i] == 0){
      data.after.base$int.time[i] <- ifelse(data.after.base$svy_interval[i] == 6, 2, data.after.base$int.time[i])
      data.after.base$int.time[i] <- ifelse(data.after.base$svy_interval[i] == 12, 3, data.after.base$int.time[i])
    }
    else{
      data.after.base$int.time[i] <- ifelse(data.after.base$svy_interval[i] == 3, 4, data.after.base$int.time[i])
      data.after.base$int.time[i] <- ifelse(data.after.base$svy_interval[i] == 6, 5, data.after.base$int.time[i])
      data.after.base$int.time[i] <- ifelse(data.after.base$svy_interval[i] == 12, 6, data.after.base$int.time[i])
    }
  }
  
  # Make intercept a factor variable
  data.after.base$int.time <- as.factor(data.after.base$int.time)
  
  # Create age/100 variable for numerical stability for control and intervention group
  data.after.base$age100 <- data.after.base$age_sample / 100
  data.after.base$age100.C <- ifelse(data.after.base[, group] == 0, 
                                         data.after.base$age_sample/ 100, 0)
  data.after.base$age100.I <- ifelse(data.after.base[, group] == 1, 
                                         data.after.base$age_sample / 100, 0)
  
  # Create indicator variable for obese BMI for control and intervention group
  data.after.base$BMI_obese <- ifelse(data.after.base$bmi_sample >= 30, 1, 0)
  data.after.base$BMI_obese.C <- ifelse(data.after.base[, group] == 0, 
                                            data.after.base$BMI_obese, 0)
  data.after.base$BMI_obese.C <- as.factor(data.after.base$BMI_obese.C)
  data.after.base$BMI_obese.I <- ifelse(data.after.base[, group] == 1, 
                                            data.after.base$BMI_obese, 0)
  data.after.base$BMI_obese.I <- as.factor(data.after.base$BMI_obese.I)
  
  # Specify variable types for hormone therapy for control and intervention group
  data.after.base$ht_sample.C <- ifelse(data.after.base[, group] == 0, 
                                            data.after.base$ht_sample, 0)
  data.after.base$ht_sample.C <- as.factor(data.after.base$ht_sample.C)
  data.after.base$ht_sample.I <- ifelse(data.after.base[, group] == 1, 
                                            data.after.base$ht_sample, 0)
  data.after.base$ht_sample.I <- as.factor(data.after.base$ht_sample.I)
  data.after.base$ht_sample <- as.factor(data.after.base$ht_sample)
  
  # Specify variable for breast cancer stage for control and intervention group
  data.after.base$breast_sample.C <- ifelse(data.after.base[, group] == 0, 
                                                data.after.base$breast_sample, 0)
  data.after.base$breast_sample.C <- as.factor(data.after.base$breast_sample.C)
  data.after.base$breast_sample.I <- ifelse(data.after.base[, group] == 1, 
                                                data.after.base$breast_sample, 0)
  data.after.base$breast_sample.I <- as.factor(data.after.base$breast_sample.I)
  data.after.base$breast_sample <- as.factor(data.after.base$breast_sample)
  
  # Change levels of survey interval to a factor with consecutive levels
  data.after.base$svy_interval <- as.factor(data.after.base$svy_interval)
  levels(data.after.base$svy_interval) <- c(1, 2, 3)
  
  # Return updated data
  return(data.after.base)
}
