#################################################################################################################################################
# This code structures the draws according to the IRT model
#################################################################################################################################################


############################################################################################################################
# Define the total number of iterations 
############################################################################################################################
n.iters.total <- nrow(draws.to.analyze.as.matrix)


############################################################################################################################
# Define the number of latent variables
############################################################################################################################
M <- 1

############################################################################################################################
# Define the structure of the model in terms of a parameters
############################################################################################################################
a.fixed.and.estimated <- matrix(NA, nrow=J, ncol=M)





############################################################################################################################
# Define empty arrays to be filled in
############################################################################################################################
a.iterations <- array(NA, c(n.iters.total, J, M))
d.iterations <- array(NA, c(n.iters.total, J))
c.iterations <- array(NA, c(n.iters.total, J))
theta.iterations <- array(NA, c(n.iters.total, n, M))


    ############################################################################################################################
    # Define 
    #   'a.fixed.and.free' 
    #   'd.fixed.and.free' 
    #   'c.fixed.and.free' 
    #   'theta.fixed.and.free'
    
    # with
    # numeric values for fixed values
    # NA's for unknown parameters
    ############################################################################################################################
    a.fixed.and.free <- matrix(NA, nrow=J, ncol=M)
    a.fixed.and.free <- as.matrix(a.fixed.and.estimated)
    d.fixed.and.free <- rep(NA, J)
    c.fixed.and.free <- rep(NA, J)

    theta.fixed.and.free <- matrix(NA, nrow=n, ncol=M)



   ############################################################################################################################
    # Loop over iterations and fill in fixed and free structure
    ############################################################################################################################
    # which.iter=1
    for(which.iter in 1:n.iters.total){

      ############################################################################################################################
      # Define the parameters for this matrix based on the fixed and free versions
      ############################################################################################################################

      a.iterations[which.iter, , ] <- a.fixed.and.free
      d.iterations[which.iter, ] <- d.fixed.and.free
      c.iterations[which.iter, ] <- c.fixed.and.free
      theta.iterations[which.iter, , ] <- theta.fixed.and.free

    } # closes loop over iterations

  ############################################################################################################################
  # Code for a 1D model
  ############################################################################################################################
  if(M==1) {
   
    ############################################################################################################################
    # Loop over iterations
    ############################################################################################################################
    # which.iter=1
    # temp <- which.iter
    # for(which.iter in temp:n.iters.total){
    if(1==0){
        for(which.iter in 1:n.iters.total){
  
        ############################################################################################################################
        # Fill in the discriminations 
        ############################################################################################################################
        for(j in 1:J){
          for(m in 1:M){    
            if(is.na(a.fixed.and.free[j,m])){
              a.iterations[which.iter,j,m] = draws.to.analyze.as.matrix[which.iter, paste("a[", j, "]", sep="")]
            }
          } # closes loop over m
        } # closes loop over j
    
        
        ############################################################################################################################
        # Fill in the locations 
        ############################################################################################################################
        for(j in 1:J){
          if(is.na(d.fixed.and.free[j])){
            d.iterations[which.iter,j] = draws.to.analyze.as.matrix[which.iter, paste("d[", j, "]", sep="")]
          }
        } # closes loop over j
      
        ############################################################################################################################
        # Fill in the lower asymptotes
        ############################################################################################################################
        for(j in 1:J){
          if(is.na(c.fixed.and.free[j])){
            c.iterations[which.iter,j] = draws.to.analyze.as.matrix[which.iter, paste("c[", j, "]", sep="")]
          }
        } # closes loop over j
        
  
        
        ############################################################################################################################
        # Fill in the latent variable values
        ############################################################################################################################
        for(i in 1:n){
          for(m in 1:M){    
            if(is.na(theta.fixed.and.free[i,m])){
              theta.iterations[which.iter,i,m] = draws.to.analyze.as.matrix[which.iter, paste("theta[", i, "]", sep="")]
            }
          } # closes loop over m
        } # closes loop over i
          
        
      
      } # closes loop over iterations

    }
     
    
    a.start=1
    a.finish=J
    colnames(draws.to.analyze.as.matrix)[a.start:a.finish]
    a.iterations[,,1] <- draws.to.analyze.as.matrix[, a.start:a.finish]
    #colMeans(a.iterations[,,1])
    
    c.start=J+1
    c.finish=J+J
    colnames(draws.to.analyze.as.matrix)[c.start:c.finish]
    c.iterations[,] <- draws.to.analyze.as.matrix[, c.start:c.finish]
    #colMeans(c.iterations[,])
    
    d.start=2*J+1
    d.finish=2*J+J
    colnames(draws.to.analyze.as.matrix)[d.start:d.finish]
    d.iterations[,] <- draws.to.analyze.as.matrix[, d.start:d.finish]
    #colMeans(d.iterations[,])
    
    theta.start=3*J+1
    theta.finish=3*J+n
    colnames(draws.to.analyze.as.matrix)[theta.start:theta.finish]
    theta.iterations[,,1] <- draws.to.analyze.as.matrix[, theta.start:theta.finish]
    #colMeans(theta.iterations[,,1])
    
    
    
  } # closes if M=1


