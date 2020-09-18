#################################################################################################################################################
# This code structures the draws according to the LCA model
#################################################################################################################################################



############################################################################################################################
# Define the number of Latent Classes
############################################################################################################################
C <- 2


############################################################################################################################
# Define empty arrays to be filled in
############################################################################################################################
gamma.iterations <- array(NA, c(n.iters.total, C))
pi.iterations <- array(NA, c(n.iters.total, C, J))
theta.iterations <- array(NA, c(n.iters.total, n))


    ############################################################################################################################
    # Define 
    #   'gamma.fixed.and.free' 
    #   'pi.fixed.and.free' 
    #   'theta.fixed.and.free'
    
    # with
    # numeric values for fixed values
    # NA's for unknown parameters
    ############################################################################################################################
    gamma.fixed.and.free <- rep(NA, C) 
    pi.fixed.and.free <- matrix(NA, nrow=C, ncol=J)
    theta.fixed.and.free <- rep(NA, n)



   ############################################################################################################################
    # Loop over iterations and fill in fixed and free structure
    ############################################################################################################################
    # which.iter=1
    for(which.iter in 1:n.iters.total){

      ############################################################################################################################
      # Define the parameters for this matrix based on the fixed and free versions
      ############################################################################################################################

      gamma.iterations[which.iter, ] <- gamma.fixed.and.free
      pi.iterations[which.iter, ,] <- pi.fixed.and.free
      theta.iterations[which.iter, ] <- theta.fixed.and.free

    } # closes loop over iterations

  
    ############################################################################################################################
    # Start plugging in iterations
    # easy because nothing really fixed
    ############################################################################################################################
    
    previous.finish = 0
    gamma.start=previous.finish+1
    gamma.finish=previous.finish+C
    colnames(draws.to.analyze.as.matrix)[gamma.start:gamma.finish]
    gamma.iterations[,] <- draws.to.analyze.as.matrix[, gamma.start:gamma.finish]
    #colMeans(a.iterations[,,1])
    
    dim(pi.iterations)
    previous.finish = gamma.finish
    pi.start=previous.finish+1
    pi.finish=previous.finish+J*C
    colnames(draws.to.analyze.as.matrix)[pi.start:pi.finish]
    for(c in 1:C){

        start.position <- pi.start + (c-1)*J
        end.position <- start.position + J-1
        colnames(draws.to.analyze.as.matrix)[start.position:end.position]
        pi.iterations[,c, ] <- draws.to.analyze.as.matrix[, start.position:end.position]  

    }

    dim(theta.iterations)
    previous.finish = pi.finish
    theta.start=previous.finish+1
    theta.finish=previous.finish+n
    colnames(draws.to.analyze.as.matrix)[theta.start:theta.finish]
    theta.iterations[,] <- draws.to.analyze.as.matrix[, theta.start:theta.finish]
    
    
    
  