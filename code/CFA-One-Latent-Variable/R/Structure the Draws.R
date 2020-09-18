#################################################################################################################################################
# This code structures the draws according to the model with one latent variable (factor)
#################################################################################################################################################


############################################################################################################################
# Define the number of latent variables
############################################################################################################################
M <- 1


############################################################################################################################
# Define the total number of iterations 
############################################################################################################################
n.iters.total <- nrow(draws.to.analyze.as.matrix)

############################################################################################################################
# Define empty arrays to be filled in
############################################################################################################################
lambda.iterations <- array(NA, c(n.iters.total, J, M))
tau.iterations <- array(NA, c(n.iters.total, J))
psi.iterations <- array(NA, c(n.iters.total, J, J))
kappa.iterations <- array(NA, c(n.iters.total, M))
phi.iterations <- array(NA, c(n.iters.total, M, M))
ksi.iterations <- array(NA, c(n.iters.total, n, M))



  ############################################################################################################################
  # Code for a one-factor model
  ############################################################################################################################
  if(M==1) {

    ############################################################################################################################
    # Define 
    #   'lambda.fixed.and.free' 
    #   'tau.fixed.and.free' 
    #   'psi.fixed.and.free' 
    #   'kappa.fixed.and.free'
    #   'phi.fixed.and.free' 
    #   'ksi.fixed.and.free' 
    
    # with
    # numeric values for fixed values
    # NA's for unknown parameters
    ############################################################################################################################
    lambda.fixed.and.free <- matrix(NA, nrow=J, ncol=M)
    lambda.fixed.and.free[1,1] = 1

    tau.fixed.and.free <- rep(NA, J)
    
    psi.fixed.and.free <- diag(NA, nrow=J, ncol=J)
    
    kappa.fixed.and.free <- rep(0, M)
    
    phi.fixed.and.free <- matrix(NA, nrow=M, ncol=M)

    ksi.fixed.and.free <- matrix(NA, nrow=n, ncol=M)
    
    
 
    
    ############################################################################################################################
    # Loop over iterations
    ############################################################################################################################
    # which.iter=1
    for(which.iter in 1:n.iters.total){

      ############################################################################################################################
      # Define the parameters for this matrix based on the fixed and free versions
      ############################################################################################################################

      lambda.iterations[which.iter, , ] <- lambda.fixed.and.free
      tau.iterations[which.iter, ] <- tau.fixed.and.free
      psi.iterations[which.iter, , ] <- psi.fixed.and.free
      kappa.iterations[which.iter, ] <- kappa.fixed.and.free
      phi.iterations[which.iter, , ] <- phi.fixed.and.free
      ksi.iterations[which.iter, , ] <- ksi.fixed.and.free
      

      ############################################################################################################################
      # Fill in the loadings 
      ############################################################################################################################
      for(j in 1:J){
        for(m in 1:M){    
          if(is.na(lambda.fixed.and.free[j,m])){
            lambda.iterations[which.iter,j,m] = draws.to.analyze.as.matrix[which.iter, paste("lambda[", j, "]", sep="")]
          }
        } # closes loop over m
      } # closes loop over j
  
      
      ############################################################################################################################
      # Fill in the intercepts 
      ############################################################################################################################
      for(j in 1:J){
        if(is.na(tau.fixed.and.free[j])){
          tau.iterations[which.iter,j] = draws.to.analyze.as.matrix[which.iter, paste("tau[", j, "]", sep="")]
        }
      } # closes loop over j
    
      
      ############################################################################################################################
      # Fill in the error variances 
      ############################################################################################################################
      for(j in 1:J){
        #for(jj in 1:J){    
          if(is.na(psi.fixed.and.free[j,j])){
            psi.iterations[which.iter,j,j] = draws.to.analyze.as.matrix[which.iter, paste("psi[", j, "]", sep="")]
          }
        #} # closes loop over jj
      } # closes loop over j

      
      ############################################################################################################################
      # Fill in the latent variable means
      ############################################################################################################################
      for(j in 1:M){
        if(is.na(kappa.fixed.and.free[m])){
          kappa.iterations[which.iter,m] = draws.to.analyze.as.matrix[which.iter, paste("kappa", sep="")]
        }
      } # closes loop over j

      
      ############################################################################################################################
      # Fill in the covariance matrix of latent variables
      ############################################################################################################################
      for(m in 1:M){    
        for(mm in 1:M){    
          if(is.na(phi.fixed.and.free[m,mm])){
            phi.iterations[which.iter,m,mm] = draws.to.analyze.as.matrix[which.iter, paste("phi", sep="")]
          }
        } # closes loop over mm
      } # closes loop over m
 
      
      ############################################################################################################################
      # Fill in the latent variable values
      ############################################################################################################################
      for(i in 1:n){
        for(m in 1:M){    
          if(is.na(ksi.fixed.and.free[i,m])){
            ksi.iterations[which.iter,i,m] = draws.to.analyze.as.matrix[which.iter, paste("ksi[", i, "]", sep="")]
          }
        } # closes loop over m
      } # closes loop over i
        
      
    
    } # closes loop over iterations
  } # closes if M=1
    

