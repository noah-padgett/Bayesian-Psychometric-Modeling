#################################################################################################################################################
# R code to compute the PPMC analyses for the Mixed-Number Subtraction BN model
#################################################################################################################################################



############################################################################################################################
# Setup a subfolder for posterior predictive datasets (and expected value matrices)
############################################################################################################################
PP.datasets.folder <- paste(PPMC.folder, "PP Datasets\\", sep="")
dir.create(PP.datasets.folder)


#################################################################################################################################################

#      		GENERATE POSTERIOR PREDICTIVE DATASETS
#         AND EXPECTED VALUE MATRICES

#         ONLY NEED TO DO THIS ONCE
#         THEN CAN READ IN SAVED OUTPUT

#################################################################################################################################################
if(1==1){
  

############################################################################################################################
# Set up arrays for 
#   posterior predicted values
#   expected values
#   useful for saving and loading workspaces
############################################################################################################################
x.postpred.array <- array(NA, c(n.iters.PPMC, n, J))
expected.values.postpred.array <- array(NA, c(n.iters.PPMC, n, J))


############################################################################################################################
# Loop over iterations
############################################################################################################################
for(which.iter in 1:n.iters.PPMC){

  print(noquote(paste("Generating Expected and Posterior Predicted Matrices for PPMC using Conditional Sums for iteration ", which.iter, sep="")))
  
  
  #################################################################################################################################################
  # Construct a pi matrix with J rows and 2 columns
  # Each row contains the pi values:
  #   the first column in pi_0
  #   the second column in pi_1
  #################################################################################################################################################
  
  pi.matrix <- matrix(NA, nrow=J, ncol=2) # Will use the +1 labeling of delta
  for(j in 1:J){
    pi.matrix[j,] <- c(
      draws.to.analyze.as.matrix[which.iter, paste("pi_0[", j, "]", sep="")], 
      draws.to.analyze.as.matrix[which.iter, paste("pi_1[", j, "]", sep="")]
    )
  }
  
  
  
  #################################################################################################################################################
  # Calculate the model-implied expected values for each data point
  #################################################################################################################################################
  expected.values <- matrix(NA, nrow=n, ncol=J)
  
  
  for(i in 1:n){
    for(j in 1:J){
      if(sum(j ==  c(4,seq(6,12),seq(14,20)))){
        expected.values[i, j] <- pi.matrix[
          j, 
          (1 + draws.to.analyze.as.matrix[which.iter, paste("delta[", i, ",", j, "]", sep="")])
        ]
      } # closes if it's an item included 
    } # closes loop over j
  } # closes loop over i
  
  

  #################################################################################################################################################
  # Construct the posterior predicted dataset
  #################################################################################################################################################
  x.postpred <- matrix(NA, nrow=n, ncol=J)
  random.number.for.post.pred <- matrix(runif(n*j), nrow=n)
  
  for(i in 1:n){
    for(j in 1:J){
      
      if(!is.na(expected.values[i,j])){
        if(expected.values[i,j] >= random.number.for.post.pred[i,j]) {x.postpred[i,j] = 1}
        else {x.postpred[i,j] = 0}
      } # closes if expected.value is not NA
    
    } # closes loop over j
  } # closes loop over i


  #################################################################################################################################################
  # Write out the 
  #   posterior predicted data
  #   model-implied expected values
  #   model-implied covariance matrix
  #   model-implied mean vector
  #################################################################################################################################################
  write.table(x.postpred, file=paste(PP.datasets.folder, "x.postpred.", which.iter, ".dat", sep=""), row.names=FALSE, col.names=FALSE)	
	write.table(expected.values, file=paste(PP.datasets.folder, "expected.values.", which.iter, ".dat", sep=""), row.names=FALSE, col.names=FALSE)	

  
  #################################################################################################################################################
  # Store the 
  #   posterior predicted data
  #   model-implied expected values
  #################################################################################################################################################
  x.postpred.array[which.iter,,] <- x.postpred
  expected.values.postpred.array[which.iter,,] <- expected.values
  
    
} # closes loop over iterations
  
  

########################################################################################################################
# Save a workspace with
#   posterior predicted data
#   expected values
########################################################################################################################
file.name <- "Posterior Predicted Data and Expected Values.Rdata"
  save(
      x.postpred.array, 
      expected.values.postpred.array,
      n.iters.PPMC,
  file=paste(PP.datasets.folder, file.name, sep="")
)



} # closes switch for generating posterior predictive datatsets and expected value matrices




#################################################################################################################################################

#  				ANALYSIS OF CONDITIONAL DISTRIBUTION OF A SUM OF A SUBSET GIVEN SUM OF OTHER SUBSET

#################################################################################################################################################

#################################################################################################################################################
# Function to compute conditional distribution table given a conditioned value of the sum of other variables
#################################################################################################################################################
conditional.distribution.sum.function <- function(
  X, 
  conditioning.variables, 
  conditioned.value, 
  outcome.variables){
  
  conditioning.sum <- rowSums(X[,conditioning.variables])
  outcome.sum <- rowSums(X[,outcome.variables])
  temp <- subset(cbind(conditioning.sum, outcome.sum), conditioning.sum==conditioned.value)
  possible.scores <- seq(0, length(outcome.variables))
  
  summary.table <- cbind(possible.scores, rep(NA, length(possible.scores)))
  colnames(summary.table) <- c("outcome.score", "count")
  for(which.score in 1:nrow(summary.table)){
    summary.table[which.score, 2] <- sum(temp[,2]==summary.table[which.score,1])
  }
  summary.table
}




#################################################################################################################################################
# Read in the
#   posterior predicted data
#   expected values
#################################################################################################################################################
file.name <- "Posterior Predicted Data and Expected Values.Rdata"
load(file=paste(PP.datasets.folder, file.name, sep=""))
  


#################################################################################################################################################
# CONDUCT PPMC OF CONDITIONAL DISTRIBUTION FOR
#   (Y_134 given Y_13)
#   ITEMS MEASURING SKILLS (1,3,4) GIVEN ITEMS MEASURING SKILLS (1,3)
#################################################################################################################################################

  #################################################################################################################################################
  # Define the 
  #   variables in the conditioning set
  #   variables in the outcome set
  #   possible scores of the latter
  #################################################################################################################################################
  conditioning.variables <- c(9, 14, 16)
  outcome.variables <- c(4, 11, 17, 18, 20)
  possible.scores <- seq(0, length(outcome.variables))
  
  
  #################################################################################################################################################
  # Loop over the 
  #   value of the conditioned sum
  #################################################################################################################################################
  # conditioned.value=0
  for(conditioned.value in 0:(length(conditioning.variables))){
    
    #################################################################################################################################################
    # Compute the realized values
    #################################################################################################################################################
    realized.summary.table <- conditional.distribution.sum.function(
      X=x, 
      conditioning.variables=conditioning.variables, 
      conditioned.value=conditioned.value, 
      outcome.variables=outcome.variables
    )
    
    #################################################################################################################################################
    # Compute the posterior predicted values
    #################################################################################################################################################
    posterior.predicted.summary.tables <- array(NA, c(n.iters.PPMC,dim(realized.summary.table)))
    
    ############################################################################################################################
    # Loop over iterations
    ############################################################################################################################
    # which.iter=1
    for(which.iter in 1:n.iters.PPMC){
      
      print(noquote(paste("Conducting PPMC using Conditional Sums for conditioned value ", conditioned.value, " of ", length(conditioning.variables), ", iteration ", which.iter, sep="")))
      
      ############################################################################################################################
      # Extract the posterior predicted dataset
      ############################################################################################################################
      x.postpred <- x.postpred.array[which.iter,,]
      
      ############################################################################################################################
      # Calculate the posterior predicted
      #   summary table
      ############################################################################################################################
      posterior.predicted.summary.tables[which.iter, , ] <- conditional.distribution.sum.function(
        X=x.postpred, 
        conditioning.variables=conditioning.variables, 
        conditioned.value=conditioned.value, 
        outcome.variables=outcome.variables
      )  
      
    } # closes loop over iterations
    
    
    #################################################################################################################################################
    # PLOT RESULTS FOR THE CONDITIONAL SCORES
    #################################################################################################################################################
    print(noquote(paste("Plotting Results from PPMC using Conditional Sums", sep="")))
    
    
    for(which.outcome.value in 1:(dim(posterior.predicted.summary.tables)[2])){
      
      #############################################################################################################################
      # Plot posterior predicted values and realized value
      # First as density
      # Then as histogram (can't get histograms to print to emf files)
      #############################################################################################################################
      values.to.plot <- posterior.predicted.summary.tables[ ,which.outcome.value, 2]
      line.to.add <- realized.summary.table[which.outcome.value,2]

      
      

      ###########################################################################
      # Define features of the plot
      ###########################################################################
      plot.width=2
      plot.height=2
     
      
      plot.margins=c(1.5,.1,.1,.1)
      
      ###########################################################################
      # Define the main name of the output files (plots of various file types)
      ###########################################################################
      file.main.name <- paste("Y_134=", realized.summary.table[which.outcome.value,1], " given Y_13=", conditioned.value, sep="")

      
            
      ###########################################################################
      # Open a window, 
      # Actually needed for creating pdf
      ###########################################################################
      if(1==1){
      windows(width=plot.width, height=plot.height)
      par(family="serif") # font
      par(mar=plot.margins)
      }
                
      plot(
          x=hist(
            values.to.plot,
            breaks = seq(
              min(values.to.plot, line.to.add)-.5,
              max(values.to.plot, line.to.add)+.5
            )
          ),
          axes=FALSE,
          ann=FALSE,
          xlim=c(
            min(
              floor(min(density(values.to.plot)$x)),
              line.to.add
            ),
            max(
              ceiling(max(density(values.to.plot)$x)),
              line.to.add
            )
          )
  

      )
      
      #############################################################################################################################
      # add line at realized value
      #############################################################################################################################
      abline(v=realized.summary.table[which.outcome.value,2], lwd=2)
      
      
      #############################################################################################################################
      # add axes
      #############################################################################################################################
      axis(side=1, font=6, labels=NA)
      axis(side = 1, lwd = 0, line = -.4, font=6)
      
 
      #############################################################################################################################
      # Copy to a .pdf
      #############################################################################################################################
      file.name <- paste(file.main.name, ".pdf", sep="")
      dev.copy2pdf(file=file.name)
      
      
    dev.off()  # closes the device opened in R
      
     
    } # closes loop over which.outcome.value   
  
  } # closes loop over conditioned.value
  

  















#################################################################################################################################################

#    			ANALYSIS OF PERSON FIT BASED ON SQUARED PEARSON RESIDUALS

#################################################################################################################################################

#################################################################################################################################################
# Read in the
#   posterior predicted data
#   expected values
#################################################################################################################################################
file.name <- "Posterior Predicted Data and Expected Values.Rdata"
load(file=paste(PP.datasets.folder, file.name, sep=""))
  


#################################################################################################################################################
# Set up arrays to capture 
# Squared Pearson residuals from realized dataset
# Person fit discrepancy measure from realized dataset

# Squared Pearson residuals from posterior predicted datasets
# Person fit discrepancy measure from posterior predicted datasets

# p.values for Squared Peason residuals
# p.values for person fit discrepancy measure

#################################################################################################################################################
realized.squared.Pearson.resids <- array(NA, c(n.iters.PPMC, n, J))
realized.person.fit <- array(NA, c(n.iters.PPMC, n))

postpred.squared.Pearson.resids <- array(NA, c(n.iters.PPMC, n, J))
postpred.person.fit <- array(NA, c(n.iters.PPMC, n))

postpred.squared.Pearson.resids.p <- array(NA, c(n, J))
postpred.person.fit.p <- array(NA, c(n))


############################################################################################################################
# Loop over iterations
############################################################################################################################
for(which.iter in 1:n.iters.PPMC){

  print(noquote(paste("Conducting PPMC using Squared Pearson and Pearson Fit Discrepancy Measurs for iteration ", which.iter, sep="")))
  
  ############################################################################################################################
  # Extract the expected value dataset
  ############################################################################################################################
  expected.values <- expected.values.postpred.array[which.iter,,]
  
  ############################################################################################################################
  # Extract in the posterior predicted dataset
  ############################################################################################################################
  x.postpred <- x.postpred.array[which.iter,,]
  
  ############################################################################################################################
  # Calculate the realized
  #   squared Pearson residuals
  #   person fit discrepancy measure
  ############################################################################################################################
  #
  realized.squared.Pearson.resids[which.iter, , ] <- as.matrix(((x - expected.values)^2)/(expected.values*(1-expected.values)))
  realized.person.fit[which.iter, ] <- sqrt(rowMeans(realized.squared.Pearson.resids[which.iter, , ], na.rm=TRUE))

  ############################################################################################################################
  # Calculate the posterior predicted
  #   squared Pearson residuals
  #   person fit discrepancy measure
  ############################################################################################################################
  #
  postpred.squared.Pearson.resids[which.iter, , ] <- as.matrix(((x.postpred - expected.values)^2)/(expected.values*(1-expected.values)))
  postpred.person.fit[which.iter, ] <- sqrt(rowMeans(postpred.squared.Pearson.resids[which.iter, , ], na.rm=TRUE))
  
  
} # closes loop over iterations
  

#############################################################################################################################
# Calculate the p-values
#############################################################################################################################
postpred.squared.Pearson.resids.p <- apply(postpred.squared.Pearson.resids > realized.squared.Pearson.resids, c(2,3), mean)
postpred.person.fit.p <- apply(postpred.person.fit > realized.person.fit, 2, mean)



#############################################################################################################################
# Write out the p-values
#############################################################################################################################

file.name <- "Person.fit.ppp.values.out"
write.table(
  cbind(seq(1,n), postpred.person.fit.p), 
  file=paste(PPMC.folder, file.name, sep=""), 
  row.names=FALSE, 
  col.names=c("examinee", "ppp.value"),
  quote=FALSE
)

