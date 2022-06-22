

generateJAGS <- function(type){
#################################################################################
# Define JAGS model for 2 variables
#################################################################################
type = as.character(type)

if(type == "2var"){

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


}

if(type == "bsplines"){
  #################################################################################
  # Define JAGS model for Bsplines
  #################################################################################
  ## nc+ni: number of covariates
  ## k+1: number of control points per spline
  ## D: degrees of the spline
  ## m+1: number of knots per spline
  ## Y: fatigue intensity
  ## X: covariates
  ## B.z1_C: basis functions for MVPA of control group
  ## B.z1_I: basis functions for MVPA of treatment group
  ## B.z2_C: basis functions for ABD of control group
  ## B.z2_I: basis functions for ABD of treatment group

  model_string <- "model {


  for (i in 1:n){
   # Likelihood (deign matrix %*% matrix of parameters to estimate)
  mu[i] <- X[i, ] %*% b + B.z1_C[i, ] %*% g1 + B.z1_I[i, ] %*% g2 + B.z2_C[i, ] %*% g3 + B.z2_I[i, ] %*% g4  ## expected response

  ## response model
    y[i] ~ dnorm(mu[i],tau)

  ## drawing the latent variables
    z1_C[i] ~ dnorm(0, taue1)
    z1_I[i] ~ dnorm(0, taue2)
    z2_C[i] ~ dnorm(0, taue3)
    z2_I[i] ~ dnorm(0, taue4)
  }

 ## creating the basis functions for the Bsplines
    for(i in 1:(M+1)){
    # knots for S(MVPA) model for control group
    knot.1_C[i] <- z1_C[knot.position[i]]

    # knots for S(MVPA) model for treatment group
    knot.1_I[i] <- z1_I[knot.position[i]]

    # knots for S(ABD) model for control group
    knot.2_C[i] <- z2_C[knot.position[i]]

    # knots for S(ABD) model for treatment group
    knot.2_I[i] <- z2_I[knot.position[i]]

  }

    for(p in 1:n){
      for(m in 1:M){
       #declaring the first layer of basis functions
        N1_C[m,1,p] <- ifelse(z1_C[p] >= knot.1_C[m] && z1_C[p] <= knot.1_C[m+1], 1, 0)
        N1_I[m,1,p] <- ifelse(z1_I[p] >= knot.1_I[m] && z1_I[p] <= knot.1_I[m+1], 1, 0)
        N2_C[m,1,p] <- ifelse(z2_C[p] >= knot.2_C[m] && z2_C[p] <= knot.2_C[m+1], 1, 0)
        N2_I[m,1,p] <- ifelse(z2_I[p] >= knot.2_I[m] && z2_I[p] <= knot.2_I[m+1], 1, 0)
      }

    #declaring the other layers of basis functions
      for(d in 2:(D+1)){
        for(m in 1:(M-d+1)){

      #S(MVPA) model for control group
        temp1[m,d,p] <-  (((z1_C[p] - knot.1_C[m]+0.0001)/(knot.1_C[m+d-1] - knot.1_C[m]+0.0003))*N1_C[m,d-1,p]
                  + ((knot.1_C[m+d] - z1_C[p]+0.0002)/(knot.1_C[m+d] - knot.1_C[m+1]+0.0004))*N1_C[m+1,d-1,p])
        N1_C[m,d,p] <- ifelse(knot.1_C[m+d-1] == knot.1_C[m],0,temp1[m,d,p])

      #S(MVPA) model for treatment group
        temp2[m,d,p] <- (((z1_I[p] - knot.1_I[m]+0.0001)/(knot.1_I[m+d-1] - knot.1_I[m]+0.0003))*N1_I[m,d-1,p]
                 + ((knot.1_I[m+d] - z1_I[p]+0.0002)/(knot.1_I[m+d] - knot.1_I[m+1]+0.0004))*N1_I[m+1,d-1,p])
        N1_I[m,d,p] <- ifelse(knot.1_C[m+d-1] == knot.1_C[m],0,temp2[m,d,p])

      #S(ABD) model for control group
        temp3[m,d,p] <- (((z2_C[p] - knot.2_C[m]+0.0001)/(knot.2_C[m+d-1] - knot.2_C[m]+0.0003))*N2_C[m,d-1,p]
                 + ((knot.2_C[m+d] - z2_C[p]+0.0002)/(knot.2_C[m+d] - knot.2_C[m+1]+0.0004))*N2_C[m+1,d-1,p])
        N2_C[m,d,p] <- ifelse(knot.1_C[m+d-1] == knot.1_C[m],0,temp3[m,d,p])

      #S(ABD) model for treatment group
        temp4[m,d,p] <- (((z2_I[p] - knot.2_I[m]+0.0001)/(knot.2_I[m+d-1] - knot.2_I[m]+0.0003))*N2_I[m,d-1,p]
                 + ( (knot.2_I[m+d] - z2_I[p]+0.0002)/(knot.2_I[m+d] - knot.2_I[m+1]+0.0004))*N2_I[m+1,d-1,p])
        N2_I[m,d,p] <- ifelse(knot.1_C[m+d-1] == knot.1_C[m],0,temp4[m,d,p])
       }

   }

    for(r in 1:(k+1)){

    #S(MVPA) model for control group
      B.z1_C[p,r] <- N1_C[r, D+1,p]

    #S(MVPA) model for treatment group
      B.z1_I[p,r] <- N1_I[r, D+1,p]

    #S(ABD) model for control group
      B.z2_C[p,r] <- N2_C[r, D+1,p]

    #S(ABD) model for treatment group
      B.z2_I[p,r] <- N2_I[r, D+1,p]
    }

  }


 ## classical measurement error model
 for(i in 1:n){
  for(j in 1:(k+1)){
    W1_C[i,j] ~ dnorm(B.z1_C[i, j]*g1[j], tauu1)
    W1_I[i,j] ~ dnorm(B.z1_I[i, j]*g2[j], tauu2)
    W2_C[i,j] ~ dnorm(B.z2_C[i, j]*g3[j] , tauu3)
    W2_I[i,j] ~ dnorm(B.z2_I[i, j]*g4[j], tauu4)

  }
 }

  ## Prior for covariates and intercept
  for (i in 1:nci) { b[i] ~ dnorm(0,0.0075) }

  ## Prior for spline parameters
  for(i in 1:(k+1)){
  g1[i] ~ dnorm(0,0.01)
  g2[i] ~ dnorm(0,0.01)
  g3[i] ~ dnorm(0,0.01)
  g4[i] ~ dnorm(0,0.01)
  }

  tau ~ dgamma(.05,.005) ## precision parameter prior
  scale <- 1/tau ## convert tau to standard GLM scale

  # prior distributions error model
  tauu1 ~  dgamma(0.01, 0.01)
  tauu2 ~  dgamma(0.01, 0.01)
  tauu3 ~  dgamma(0.01, 0.01)
  tauu4 ~  dgamma(0.01, 0.01)
  taue1 ~  dgamma(0.01, 0.01)
  taue2 ~  dgamma(0.01, 0.01)
  taue3 ~  dgamma(0.01, 0.01)
  taue4 ~  dgamma(0.01, 0.01)

  sigma_u1 <- 1/tauu1
  sigma_u2 <- 1/tauu2
  sigma_u3 <- 1/tauu3
  sigma_u4 <- 1/tauu4

}"
}

if(type == "penalized"){

  #################################################################################
  # Define JAGS model for penalized splines
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

  sigma_u1 <- 1/tauu1
  sigma_u2 <- 1/tauu2
  sigma_u3 <- 1/tauu3
  sigma_u4 <- 1/tauu4
}"


}
return("model_string" = model_string)
}
