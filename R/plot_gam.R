################################################################################################
# Function to plot factor by curve model 
################################################################################################
## Variables
# fit = factor by curve model with varying intercept
# pa_vars = vector of PA variable names modeled as smooth functions
# covars = vector of covariate variable names modeled as linear predictors
# time = vector of time periods and group indicators
# y = outcome variable name
# ymin and ymax = min and max values for plotting
# data = data used to fit model
# bounds = logical value to determine whether you want variability bands when plotting
# group = variable for grouping factor

plot_FBC <- function(fit, pa_vars, covars, time, y, ymin, ymax, data, bounds = F, plot.title, group){
  
  # Specify number of grid points
  ngrid <- 101
  
  # List for storing plotting variables
  pa.list <- list()
  for(pa in pa_vars){
    pa.list[[pa]] <- data.frame(matrix(NA, nrow = ngrid*6, ncol = 6))
    colnames(pa.list[[pa]]) <- c("x", "s.x", "upper", "lower", "type", "group")
  }
  
  # Index value for loop
  j <- 1
  
  # Loop over each time period
  for(t in time){
    # If time period in 1-3, group type = control
    if(t %in% 1:3){
      type <- 0
    }
    # If time period in 4-6, group type = treatment
    else{
      type <- 1
    }
    
    # Create index sequence for how many variables to predict
    ind <- ((j-1)*ngrid+1):(j*ngrid)
    
    ## Create new dataframe for plotting
    
    # Obtain factor covariates
    factor.vars <- names(data)[ sapply(data, is.factor) ]
    factor.covars <- factor.vars[which(factor.vars %in% covars)]
    
    # Obtain continuous covariates
    cont.covars <- covars[which(!(covars %in% factor.covars))]
    
    # Set covariates to their mean/median values
    cov.avg <- c(sapply(data[which(data[ , group] == type), factor.covars, drop = F], 
                        function(x) median(as.numeric(as.character(x)))), 
                 sapply(data[which(data[ , group] == type), cont.covars, drop = F], mean))
    
    # Create matrix of averaged covariates with ngrid rows
    cov.mat <- matrix(rep(cov.avg, each = ngrid), nrow = ngrid)
    
    # Take mean value of PA variables
    pa.means <- rep(colMeans(data[ , pa_vars, drop = F]), each = ngrid)
    #pa.means <- rep(colMeans(data[which(data[ , group] == type), pa_vars, drop = F]), each = ngrid)
    
    # Create matrix of averaged PA variables with ngrid rows
    pa.mat <- matrix(pa.means, ncol = length(pa_vars), nrow = ngrid)
    
    ## Combine PA variables and covairates into one dataframe
    newdata <- data.frame(cbind(pa.mat, rep(type, length = ngrid), rep(t, length = ngrid), cov.mat))
    colnames(newdata) <- c(pa_vars, group, "int.time", factor.covars, cont.covars)
    
    # Loop over each PA variable for plotting
    for(pa in pa_vars){
      
      # Create grid of PA variable for plotting
      type.ind <- which(data[ , group] == type)
      pa.grid <- seq(min(data[type.ind, pa]), 
                     max(data[type.ind, pa]), length = ngrid)
      newdata.tmp <- newdata
      newdata.tmp[ , pa] <- pa.grid
      
      # Get predicted values for pa variable + CI
      pred.pa <- predict(fit$gam, newdata = newdata.tmp, se = T)
      
      # Rescale the response to original units of measure
      s.pa <- pred.pa$fit * sd(data[ , y]) + mean(data[ , y])
      s.pa.lw <- s.pa - qnorm(0.975) * pred.pa$se.fit * sd(data[ , y])
      s.pa.up <- s.pa + qnorm(0.975) * pred.pa$se.fit * sd(data[ , y])
      pa.orig <- pa.grid * sd(data[ , unlist(strsplit(pa, "_sc"))]) +
        mean(data[ , unlist(strsplit(pa, "_sc"))])
      # s.pa <- pred.pa$fit * sd(data[which(data[ , group] == type), y]) + mean(data[which(data[ , group] == type), y])
      # s.pa.lw <- s.pa - qnorm(0.975) * pred.pa$se.fit * sd(data[which(data[ , group] == type), y])
      # s.pa.up <- s.pa + qnorm(0.975) * pred.pa$se.fit * sd(data[which(data[ , group] == type), y])
      # pa.orig <- pa.grid * sd(data[which(data[ , group] == type), unlist(strsplit(pa, "_sc"))]) +
      #   mean(data[data[ , group] == type, unlist(strsplit(pa, "_sc"))])
      
      # Combine results into dataframe
      pa.list[[pa]][ind, ] <- cbind(pa.orig, s.pa, s.pa.lw, s.pa.up, t, type)
    }
    
    # Update loop count
    j <- j + 1
  }
  
  # List for plots
  plotlist <- list()
  j <- 1 # index for loop
  
  # Loop over each PA variable
  for(pa in pa_vars){
    
    # Get unit of measure for pa variable
    xlab <- paste0("Difference in ", plot.title[j], " (min/day)")
    
    # Get unit of measure for response
    if(grepl("intensity", y)){
      ylab <- "Difference in fatigue intensity after baseline"
    }
    else{
      ylab <- "Difference in fatigue interference after baseline"
    }
    
    # Plot pa variable
    est.df <- pa.list[[pa]]
    
    # Create factor for group (control of tret) for rug plot
    est.df$group <- as.factor(est.df$group)
    
    # Create factor for time period and group type
    est.df$type <- as.factor(est.df$type)
    levels(est.df$type) <- c("Ctrl-M3", "Ctrl-M6", "Ctrl-M12", "Trt-M3", "Trt-M6", "Trt-M12")
    
    # Get original x values for plotting
    x_orig <- unlist(strsplit(pa, "_sc"))
    est.df$x_orig <- rep(data[,x_orig], length = nrow(est.df))
    
    # Plot pa variable
    plotlist[[pa]] <- ggplot(est.df, aes(x = x, y = s.x, group = type, color = type)) +
      geom_line() + theme_bw() + ylim(ymin, ymax) + 
      labs(x = xlab, y = ylab, title = plot.title[j]) +
      #geom_rug(aes(x = x_orig, color = group), sides = "b") +
      theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"), 
            text = element_text(size = 10), plot.subtitle = element_text(hjust = 0.5))
    
    # If plotting with confidence intervals
    if (bounds){
      plotlist[[pa]] <- plotlist[[pa]] + geom_line(aes(x = x, y = lower, color = type), linetype = "dashed") +
        geom_line(aes(x = x, y = upper, color = type), linetype = "dashed") 
    }
    
    # Update loop count
    j <- j + 1
  }
  
  return(plotlist)
}
