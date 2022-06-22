dataPostprocessingLin = function(samples.jags){
  ##matrix to store the g estimates
  g_mat = matrix(NA, 1000,4)
  g_mat[,1] <- samples.jags$g1[1, ,1]
  g_mat[,2] <- samples.jags$g2[1, ,1]
  g_mat[,3] <- samples.jags$g3[1, ,1]
  g_mat[,4] <- samples.jags$g4[1, ,1]
  g_mat_names <- c("MVPA_C", "MVPA_I","ABD_C","ABD_I")

  ##creating the trace plots of the g estimates and saving them
  gamma.list <- lapply(1:4, function(i) {
    ggplot(data = as.data.frame(g_mat[ ,i]), aes(x = 1:1000, y = g_mat[ ,i] ))+
      geom_line() + xlab("Index") +  theme_bw()+
      ggtitle(g_mat_names[i])+ylab(NULL)
  })
  png(file="./plots/2var/gamma.png")
  ggarrange(plotlist = gamma.list, nrow = 2, ncol = 2)
  dev.off()

  ##matrix to store the beta estimates
  beta_mat = matrix(NA, 1000, ncol(X.mat))
  for(i in 1:ncol(X.mat)){
    beta_mat[ , i] =  samples.jags$beta[i, ,1]
  }
  ##creating the trace plots of the beta estimates and saving them
  beta.list <- lapply(1:ncol(X.mat), function(i) {
    ggplot(data = as.data.frame(beta_mat[ ,i]), aes(x = 1:1000, y = beta_mat[ ,i] ))+
      geom_line() + xlab("Index") +  theme_bw()+
      ggtitle(colnames(X.mat)[i])+ylab(NULL)
  })

  png(file="./plots/2var/beta.png")
  ggarrange(plotlist = beta.list, nrow = 4, ncol = 5)
  dev.off()

  ##creating the trace plots of the sigma estimates and saving them
  sigma_mat = matrix(NA, 1000, 5)
  sigma_mat[ ,1] <- samples.jags$sigma_u1[1, ,1]
  sigma_mat[ ,2] <- samples.jags$sigma_u2[1, ,1]
  sigma_mat[ ,3] <- samples.jags$sigma_u3[1, ,1]
  sigma_mat[ ,4] <- samples.jags$sigma_u4[1, ,1]
  sigma_mat[ ,5] <- samples.jags$sigma_y[1, ,1]
  s2_names <- c("sigma_u1 for MVPA_C", "sigma_u2 for MVPA_I","sigma_u3 for ABD_C", "sigma_u4 for ABD_I", "sigma_y")

  sigmas.list <- lapply(1:5, function(i) {
    ggplot(data = as.data.frame(sigma_mat[ ,i]), aes(x = 1:1000, y = sigma_mat[ ,i] ))+
      geom_line() + xlab("Index")+  theme_bw() +
      ggtitle(s2_names[i])+ylab(NULL)
  })
  png(file="./plots/2var/sigma_plots.png")
  ggarrange(plotlist = sigmas.list, nrow = 3, ncol = 2)
  dev.off()

  return(list("gamma.list" = gamma.list,"beta.list"=beta.list, "sigmas.list"=sigmas.list ))
}
