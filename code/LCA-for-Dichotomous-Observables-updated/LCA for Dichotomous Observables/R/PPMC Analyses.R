#################################################################################################################################################
# R code to compute PPMC Analyses for LCA
#################################################################################################################################################




#################################################################################################################################################

#      		GENERATE POSTERIOR PREDICTIVE DATASETS
#         AND EXPECTED VALUE MATRICES
#################################################################################################################################################

############################################################################################################################
# Setup a subfolder for posterior predictive datasets (and expected value matrices)
############################################################################################################################
PP.datasets.folder <- paste(PPMC.folder, "PP Datasets\\", sep="")
dir.create(PP.datasets.folder)



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
# which.iter=1
for(which.iter in 1:n.iters.PPMC){

  print(noquote(paste("Generating Expected and Posterior Predicted Matrices for PPMC for iteration ", which.iter, sep="")))
  
  
  
  #################################################################################################################################################
  # Calculate the model-implied expected values for each data point
  #################################################################################################################################################
  expected.values <- matrix(NA, nrow=n, ncol=J)
  
  
  for(i in 1:n){
    for(j in 1:J){
      #if(sum(j ==  c(4,seq(6,12),seq(14,20)))){
        expected.values[i, j] <- pi.iterations.for.PPMC[which.iter, theta.iterations.for.PPMC[which.iter,i], j]      
      #} # closes if it's an item included 
    } # closes loop over j
  } # closes loop over i
  
  

  #################################################################################################################################################
  # Construct the posterior predicted dataset
  #################################################################################################################################################
  x.postpred <- (expected.values >= matrix(runif(n*j), nrow=n))*1
  

  #################################################################################################################################################
  # Write out the 
  #   posterior predicted data
  #   model-implied expected values
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



#################################################################################################################################################

#  				ANALYSIS OF PERSON AND OBSERVABLE FIT BASED ON SQUARED PEARSON RESIDUALS

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
# Observalbe fit discrepancy measure from realized dataset

# Squared Pearson residuals from posterior predicted datasets
# Person fit discrepancy measure from posterior predicted datasets
# Observalbe fit discrepancy measure from posterior predicted datasets

# p.values for Squared Peason residuals
# p.values for person fit discrepancy measure
# p.values for observable fit discrepancy measure

#################################################################################################################################################
realized.squared.Pearson.resids <- array(NA, c(n.iters.PPMC, n, J))
realized.person.fit <- array(NA, c(n.iters.PPMC, n))
realized.observable.fit <- array(NA, c(n.iters.PPMC, J))

postpred.squared.Pearson.resids <- array(NA, c(n.iters.PPMC, n, J))
postpred.person.fit <- array(NA, c(n.iters.PPMC, n))
postpred.observable.fit <- array(NA, c(n.iters.PPMC, J))

postpred.squared.Pearson.resids.p <- array(NA, c(n, J))
postpred.person.fit.p <- array(NA, c(n))
postpred.observable.fit.p <- array(NA, c(J))


############################################################################################################################
# Loop over iterations
############################################################################################################################
for(which.iter in 1:n.iters.PPMC){

  print(noquote(paste("Conducting PPMC using Squared Pearson, Pearson Fit, and Observable Fit Discrepancy Measurs for iteration ", which.iter, sep="")))
  
  ############################################################################################################################
  # Extract the expected value dataset
  ############################################################################################################################
  expected.values <- expected.values.postpred.array[which.iter,,]
  
  ############################################################################################################################
  # Extract the posterior predicted dataset
  ############################################################################################################################
  x.postpred <- x.postpred.array[which.iter,,]
    
  ############################################################################################################################
  # Calculate the realized
  #   squared Pearson residuals
  #   person fit discrepancy measure
  #   observable fit discrepancy measure
  ############################################################################################################################
  #
  realized.squared.Pearson.resids[which.iter, , ] <- as.matrix(((x - expected.values)^2)/(expected.values*(1-expected.values)))
  realized.person.fit[which.iter, ] <- sqrt(rowMeans(realized.squared.Pearson.resids[which.iter, , ], na.rm=TRUE))
  realized.observable.fit[which.iter, ] <- sqrt(colMeans(realized.squared.Pearson.resids[which.iter, , ], na.rm=TRUE))

  ############################################################################################################################
  # Calculate the posterior predicted
  #   squared Pearson residuals
  #   person fit discrepancy measure
  #   observable fit discrepancy measure
  ############################################################################################################################
  #
  postpred.squared.Pearson.resids[which.iter, , ] <- as.matrix(((x.postpred - expected.values)^2)/(expected.values*(1-expected.values)))
  postpred.person.fit[which.iter, ] <- sqrt(rowMeans(postpred.squared.Pearson.resids[which.iter, , ], na.rm=TRUE))
  postpred.observable.fit[which.iter, ] <- sqrt(colMeans(postpred.squared.Pearson.resids[which.iter, , ], na.rm=TRUE))
  
  
} # closes loop over iterations
  

#############################################################################################################################
# Calculate the p-values
#############################################################################################################################
postpred.squared.Pearson.resids.p <- apply(postpred.squared.Pearson.resids > realized.squared.Pearson.resids, c(2,3), mean)
postpred.person.fit.p <- apply(postpred.person.fit > realized.person.fit, 2, mean)
postpred.observable.fit.p <- apply(postpred.observable.fit > realized.observable.fit, 2, mean)



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

file.name <- "Observable.fit.ppp.values.out"
write.table(
  cbind(seq(1,J), postpred.observable.fit.p), 
  file=paste(PPMC.folder, file.name, sep=""), 
  row.names=FALSE, 
  col.names=c("observable", "ppp.value"),
  quote=FALSE
)


