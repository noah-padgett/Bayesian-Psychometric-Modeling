#################################################################################################################################################
# R code to compute the PPMC analyses for the CFA model
#################################################################################################################################################




#################################################################################################################################################

#      		GENERATE POSTERIOR PREDICTIVE DATASETS
#         & MODEL-IMPLIED MEAN AND COVARIANCE MATRICES

#         ONLY NEED TO RUN THIS ONCE, THEN CAN REVISIT LATER CODE
#################################################################################################################################################

############################################################################################################################
# Setup a subfolder for posterior predictive datasets
############################################################################################################################
PP.datasets.folder <- paste(PPMC.folder, "PP Datasets\\", sep="")
dir.create(PP.datasets.folder)

  
############################################################################################################################
# Loop over iterations
############################################################################################################################
for(which.iter in 1:n.iters.total){

  
  
  #################################################################################################################################################
  # Extract the values of the needed model parameters for this iteration
  #################################################################################################################################################
  lambda <- matrix(lambda.iterations[which.iter, , ], nrow=J)
  tau <- tau.iterations[which.iter, ]
  psi <- psi.iterations[which.iter, , ]
  kappa <- kappa.iterations[which.iter, ]
  phi <- matrix(phi.iterations[which.iter, , ], nrow=M)
  ksi <- matrix(ksi.iterations[which.iter, , ], nrow=n)
 

  
  #################################################################################################################################################
  # Calculate the model-implied expected values for each data point
  #################################################################################################################################################
  expected.values <- (matrix(rep(1, n), nrow=n) %*% tau) + (ksi %*% t(lambda))

  
  #################################################################################################################################################
  # Simulate error values in a matrix
  #################################################################################################################################################
  error.values <- matrix(NA, nrow=n, ncol=J)
  for(j in 1:J){
		error.values[, j] <- rnorm(n=n, mean=0, sd=sqrt(diag(psi)[j]))
	}

  
  #################################################################################################################################################
  # Construct the posterior predicted dataset
  #################################################################################################################################################
  x.postpred <- expected.values + error.values
  


  #################################################################################################################################################
  # Compute the model-implied mean vector
  #################################################################################################################################################
  model.implied.mean.vector <- tau + lambda %*% kappa


  #################################################################################################################################################
  # Compute the model-implied covariance matrix
  #################################################################################################################################################
  model.implied.covariance.matrix <- lambda %*% phi %*% t(lambda) + psi

  
  #################################################################################################################################################
  # Write out the 
  #   posterior predicted data
  #   model-implied expected values
  #   model-implied covariance matrix
  #   model-implied mean vector
  #################################################################################################################################################
  write.table(x.postpred, file=paste(PP.datasets.folder, "x.postpred.", which.iter, ".dat", sep=""), row.names=FALSE, col.names=FALSE)	
	write.table(expected.values, file=paste(PP.datasets.folder, "expected.values.", which.iter, ".dat", sep=""), row.names=FALSE, col.names=FALSE)	
	write.table(model.implied.covariance.matrix, file=paste(PP.datasets.folder, "model.implied.covariance.matrix.iteration.", which.iter, ".dat", sep=""), row.names=FALSE, col.names=FALSE)	
	write.table(model.implied.mean.vector, file=paste(PP.datasets.folder, "model.implied.mean.vector.iteration.", which.iter, ".dat", sep=""), row.names=FALSE, col.names=FALSE)	

  
} # closes loop over iterations
  
  


#################################################################################################################################################

#  				ANALYSIS OF CORRELATIONS

#################################################################################################################################################

#################################################################################################################################################
# Set up arrays to capture 

# correlation matrices from posterior predicted data sets
# p.values for correlations
#################################################################################################################################################
postpred.cors <- array(NA, c(n.iters.total, J, J))
postpred.cors.p <- array(NA, c(J, J))



############################################################################################################################
# Loop over iterations
############################################################################################################################
# which.iter=1
for(which.iter in 1:n.iters.total){

  ############################################################################################################################
  # Read in the posterior predicted dataset
  ############################################################################################################################
  x.postpred <- read.table(file=paste(PP.datasets.folder, "x.postpred.", which.iter, ".dat", sep=""))

  ############################################################################################################################
  # Calculate the posterior predicted
  #   correlations
  ############################################################################################################################
  postpred.cors[which.iter, , ] <- cor(x.postpred)
  

} # closes loop over iterations
  



#################################################################################################################################################
# PLOT RESULTS FOR THE CORRELATIONS
#################################################################################################################################################

###########################################################################
# Define features of the plot
###########################################################################
plot.width=3.5
plot.height=3.5
line.width=1

plot.margins=c(1.2, .6, .2, .6)



###########################################################################
# Define the main name of the output files (plots of various file types)
###########################################################################
file.main.name <- paste("Plot of Cors (below)", sep="")



#################################################################################################################################################
# CREATE PLOT
#################################################################################################################################################

###########################################################################
# Open a window, for creating the plot, which will then be saved as file types
###########################################################################

windows(width=plot.width, height=plot.height)
par(family="serif") # font
par(mar=plot.margins)


mat=matrix(seq(1:(J*J)), ncol=J)
layout(mat=mat, 
       widths=rep(1, ncol(mat)),
       heights=rep(1, nrow(mat)), 
       respect=FALSE
)



#############################################################################################################################
# Plot the correlations below the main diagonal
# Plot the variable labels along the diagonal
#############################################################################################################################

var.labels <- colnames(x)


for(j in 1:J){
  for(jj in 1:J){
    
    
    
    #############################################################################################################################
    # Plot the nothing above the main diagonal
    #############################################################################################################################
    if(j>jj){
      
      #############################################################################################################################
      # Plot density of posterior predicted values
      #############################################################################################################################
      values.to.plot <- postpred.cors[ ,j, jj]
      vertical.line.to.plot <- cor(x)[j,jj]
      plot(
        x=1,#density(values.to.plot),
        ann=FALSE,
        axes=FALSE,
        type="n"
      )
    }

    
    
    #############################################################################################################################
    # Plot the correlations below the main diagonal
    #############################################################################################################################
    if(j<jj){
      
      #############################################################################################################################
      # Plot density of posterior predicted values
      #############################################################################################################################
      values.to.plot <- postpred.cors[ ,j, jj] 
      vertical.line.to.plot <- cor(x)[j,jj]
      plot(
        density(values.to.plot),
        ann=FALSE,
        xlab=NA,
        ylab=NA,
        axes=F,
        xlim=c(
          min(floor(100*min(density(values.to.plot)$x))/100, vertical.line.to.plot),
          max(ceiling(100*max(density(values.to.plot)$x))/100, vertical.line.to.plot)
        ),
        ylim=c(
          0, 
          ceiling(max(density(values.to.plot)$y))
        ),
        font=6, font.axis=6, font.main=6, font.lab=6
      )
      
      axis(side=1, font=6, tcl=-.2, mgp=c(1,.2,0))
      #############################################################################################################################
      # add line at realized value
      #############################################################################################################################
      abline(v=vertical.line.to.plot)
    }    
    
    
    
    #############################################################################################################################
    # Plot the variable labels along the diagonal
    #############################################################################################################################
    if(j==jj){
      plot(x=c(0,0), y=c(0,0),
           type="n",
           ann=FALSE,
           axes=FALSE,
           xlim=c(-1, 1),
           ylim=c(-1, 1),
           font=6, font.axis=6, font.main=6, font.lab=6
      )
      text(x=0, y=0, bquote(italic(.(var.labels[j]))), font=6, cex=1.5)
    }
    
  } # closes first loop over plotting 
} # closes second loop over plotting 



#############################################################################################################################
# Copy to a .pdf
#############################################################################################################################
file.name <- paste(output.folder, file.main.name, ".pdf", sep="")
dev.copy2pdf(file=file.name)

dev.off()  # closes the device opened in R





#############################################################################################################################
# Calculate the p-values
#############################################################################################################################

for(j in 1:J){
  for(jj in j:J){
 
		#############################################################################################################################
		# Define the realized and posterior predicted values
		#############################################################################################################################
		postpred.values <- postpred.cors[ ,j,jj]
		realized.values <- cor(x)[j,jj]

		#############################################################################################################################
		# Calculate the p-value
		#############################################################################################################################	
		postpred.cors.p[j,jj] <- mean(postpred.values > realized.values)

	} # closes first loop over plotting of marginal covariances
} # closes second loop over plotting of marginal covariances

#############################################################################################################################
# Write out the p-values
#############################################################################################################################

file.name <- "Correlations.ppp.values.out"
write.table(postpred.cors.p, file=paste(PPMC.folder, file.name, sep=""), row.names=FALSE, col.names=colnames(x))




#################################################################################################################################################

#      		ANALYSIS OF LR (aka the ML fit function)

#################################################################################################################################################

#################################################################################################################################################
# Set up arrays to capture 
#  realized LR fit
#	 posterior predicted LR fit

#	p.values for LR fit
#################################################################################################################################################
realized.LR.fit <- array(NA, c(n.iters.total))
postpred.LR.fit <- array(NA, c(n.iters.total))

postpred.LR.fit.p <- NA


#################################################################################################################################################
# Define the fit function
#################################################################################################################################################

LR.function <- function(data.cov.matrix, mod.imp.cov.matrix, n){

  F.ML <- log(det(mod.imp.cov.matrix)) + sum(diag((data.cov.matrix %*% solve(mod.imp.cov.matrix)))) - log(det(data.cov.matrix)) - ncol(data.cov.matrix)
	LR <- (n-1)*F.ML
	LR

}



############################################################################################################################
# Loop over iterations
############################################################################################################################

for(which.iter in 1:n.iters.total){

  ############################################################################################################################
  # Read in the 
  #   posterior predicted dataset
  #   model-implied covariance matrix
  ############################################################################################################################
  x.postpred <- read.table(file=paste(PP.datasets.folder, "x.postpred.", which.iter, ".dat", sep=""))
  model.implied.covariance.matrix <- as.matrix(read.table(file=paste(PP.datasets.folder, "model.implied.covariance.matrix.iteration.", which.iter, ".dat", sep="")))
  
  
  ############################################################################################################################
  # Calculate the 
  #   realized LR value
  #   posterior predicted LR value
  ############################################################################################################################
  realized.LR.fit[which.iter] <- LR.function(data.cov.matrix=cov(x), mod.imp.cov.matrix=model.implied.covariance.matrix, n=n)
  postpred.LR.fit[which.iter] <- LR.function(data.cov.matrix=cov(x.postpred), mod.imp.cov.matrix=model.implied.covariance.matrix, n=n)
  
} # closes loop over iterations
  


#################################################################################################################################################
# PLOT RESULTS
#################################################################################################################################################

#############################################################################################################################
# Define the realized and posterior predicted values
#############################################################################################################################
realized.values <- realized.LR.fit
postpred.values <- postpred.LR.fit


###########################################################################
#  Define features of the plots 
###########################################################################
plot.width=3
plot.height=3
plot.margins=c(1.5,1.5,0,0)+.1
if(write.axis.labels==TRUE){
  plot.margins=plot.margins+c(1,1,0,0)
} # closes switch for if writing axis labels

point.size=.6 
point.type=19
xlimit<-range(realized.values,postpred.values)
ylimit<-xlimit


###########################################################################
#  Define main name of output file
###########################################################################
file.main.name <- paste("Scatterplot of LR", sep="")


###########################################################################
# Open a window for creating the plot, which will then be saved as file types
###########################################################################

windows(width=plot.width, height=plot.height)
par(family="serif") # font
par(mar=plot.margins)
  
plot(
  x=realized.values,
  y=postpred.values,
  xlim=xlimit,
  ylim=ylimit,
  cex=point.size,
  pch=point.type,
  axes=FALSE,
  ann=FALSE
)

#################################################################################################################################################
# Add the unit line
#################################################################################################################################################
abline(a=0, b=1)
box()
axis(side=1, font=6, labels=NA)
axis(side = 1, lwd = 0, line = -.4, font=6)
axis(side=2, font=6, labels=NA)
axis(side = 2, lwd = 0, line = -.4, font=6)

if(write.axis.labels==TRUE){
  mtext(side = 1, line = 1.5, text="Realized LR")
  mtext(side = 2, line = 1.5, text="Posterior Predicted LR")
} # closes switch for if writing axis labels

#############################################################################################################################
# Copy to a .pdf
#############################################################################################################################
file.name <- paste(output.folder, file.main.name, ".pdf", sep="")
dev.copy2pdf(file=file.name)

dev.off()  # closes the device opened in R




#############################################################################################################################
# Calculate the p-values
#############################################################################################################################
postpred.LR.fit.p <- mean(postpred.values > realized.values)


#############################################################################################################################
# Write out the p-values
#############################################################################################################################
file.name <- "LR.fit.ppp.value.out"
write.table(postpred.LR.fit.p, file=paste(PPMC.folder, file.name, sep=""), row.names=FALSE, col.names=FALSE)








#################################################################################################################################################

#          ANALYSIS OF SRMR

#################################################################################################################################################


#################################################################################################################################################
# Set up arrays to capture 
#  realized SRMR
#	 posterior predicted SRMR

#	p.values for LR fit
#################################################################################################################################################
realized.SRMR <- array(NA, c(n.iters.total))
postpred.SRMR <- array(NA, c(n.iters.total))

postpred.SRMR.p <- NA


#################################################################################################################################################
# Define the SRMR
#################################################################################################################################################

SRMR.function <- function(data.cov.matrix, mod.imp.cov.matrix){

  J=nrow(data.cov.matrix)	
	temp <- matrix(NA, nrow=J, ncol=J)
	
	for(j in 1:J){
		for(jprime in 1:j){
			temp[j, jprime] <- ((data.cov.matrix[j, jprime] - mod.imp.cov.matrix[j, jprime])/(data.cov.matrix[j, j] * data.cov.matrix[jprime, jprime]))^2
		}
	}
  
	SRMR <- sqrt((2*sum(temp,na.rm=TRUE))/(J*(J+1)))
	SRMR

}



############################################################################################################################
# Loop over iterations
############################################################################################################################
for(which.iter in 1:n.iters.total){

  ############################################################################################################################
  # Read in the 
  #   posterior predicted dataset
  #   model-implied covariance matrix
  ############################################################################################################################
  x.postpred <- read.table(file=paste(PP.datasets.folder, "x.postpred.", which.iter, ".dat", sep=""))
  model.implied.covariance.matrix <- as.matrix(read.table(file=paste(PP.datasets.folder, "model.implied.covariance.matrix.iteration.", which.iter, ".dat", sep="")))
  
  
  ############################################################################################################################
  # Calculate the 
  #   realized SRMR value
  #   posterior predicted SRMR value
  ############################################################################################################################
  realized.SRMR[which.iter] <- SRMR.function(data.cov.matrix=cov(x), mod.imp.cov.matrix=model.implied.covariance.matrix)
  postpred.SRMR[which.iter] <- SRMR.function(data.cov.matrix=cov(x.postpred), mod.imp.cov.matrix=model.implied.covariance.matrix)
  
} # closes loop over iterations
  


#################################################################################################################################################
# PLOT RESULTS
#################################################################################################################################################

  #############################################################################################################################
  # Define the realized and posterior predicted values
	#############################################################################################################################
	realized.values <- realized.SRMR
	postpred.values <- postpred.SRMR


  ###########################################################################
  #  Define main name of output file
  ###########################################################################
  file.main.name <- paste("Realized SRMR", sep="")


  ###########################################################################
  #  Define features of the plots 
  ###########################################################################
  plot.width=3
  plot.height=2.5
  plot.margins=c(1.5,.5,0,.5)+.1
  if(write.axis.labels==TRUE){
    plot.margins=plot.margins+c(1,0,0,0)
  } # closes switch for if writing axis labels
  

  ###########################################################################
  # Open a window, for creating the plot, which will then be saved as file types
  ###########################################################################
  
  windows(width=plot.width, height=plot.height)
  par(family="serif") # font
  par(mar=plot.margins)
  
  plot(
    density(realized.values),
    cex=point.size,
    axes=FALSE,
    ann=FALSE
  )
  
  #################################################################################################################################################
  # Add the axes
  #################################################################################################################################################
  axis(side=1, font=6, labels=NA)
  axis(side = 1, lwd = 0, line = -.4, font=6)
  
  
  if(write.axis.labels==TRUE){
    mtext(side = 1, line = 1.5, text="Realized SRMR")
  } # closes switch for if writing axis labels

  
  #############################################################################################################################
  # Copy to a .pdf
  #############################################################################################################################
  file.name <- paste(output.folder, file.main.name, ".pdf", sep="")
  dev.copy2pdf(file=file.name)
  
  dev.off()  # closes the device opened in R
    
  
 
	#############################################################################################################################
	# Summarize the realized values of SRMR
	#############################################################################################################################
	temp <- as.mcmc(realized.values)
  summary.stats <- summary(temp)
  probability.for.HPD=.95


  summary.statistics <- cbind(
    matrix(summary.stats$statistics, ncol=4), 
    matrix(summary.stats$quantiles, ncol=5), 
    matrix(HPDinterval(temp, prob=probability.for.HPD)[1,], ncol=2)
  )

  colnames(summary.statistics) <- c(
    row.names(as.matrix(summary.stats$statistics)),
    row.names(as.matrix(summary.stats$quantiles)),
    c("95% HPD lower", "95% HPD Upper")
  )


  #########################################################################################################################
  # Write out the summary statistics
  #########################################################################################################################
  file.name <- "Realized.SRMR.summary.csv"
  write.csv(
    x=summary.statistics,
    file=paste(PPMC.folder, file.name, sep="")
  )


#############################################################################################################################
# Define the plot name and open it 
#############################################################################################################################

###########################################################################
#  Define main name of output file
###########################################################################
file.main.name <- paste("Scatterplot of SRMR", sep="")


###########################################################################
#  Define features of the plots 
###########################################################################
plot.width=3
plot.height=3
plot.margins=c(1.5,1.5,0,0)+.1
if(write.axis.labels==TRUE){
  plot.margins=plot.margins+c(1,1,0,0)
} # closes switch for if writing axis labels

point.size=.6
point.type=19
xlimit<-range(realized.values,postpred.values)
ylimit<-xlimit


###########################################################################
# Open a window, for creating the plot, which will then be saved as file types
###########################################################################

windows(width=plot.width, height=plot.height)
par(family="serif") # font
par(mar=plot.margins)

plot(
  x=realized.values,
  y=postpred.values,
  xlim=xlimit,
  ylim=ylimit,
  cex=point.size,
  pch=point.type,
  axes=FALSE,
  ann=FALSE
)


#################################################################################################################################################
# Add the unit line
#################################################################################################################################################
abline(a=0, b=1)
box()
axis(side=1, font=6, labels=NA)
axis(side = 1, lwd = 0, line = -.4, font=6)
axis(side=2, font=6, labels=NA)
axis(side = 2, lwd = 0, line = -.4, font=6)

if(write.axis.labels==TRUE){
  mtext(side = 1, line = 1.5, text="Realized SRMR")
  mtext(side = 2, line = 1.5, text="Posterior Predicted SRMR")
} # closes switch for if writing axis labels

#############################################################################################################################
# Copy to a .pdf
#############################################################################################################################
file.name <- paste(output.folder, file.main.name, ".pdf", sep="")
dev.copy2pdf(file=file.name)

dev.off()  # closes the device opened in R

}




#############################################################################################################################
# Calculate the p-values
#############################################################################################################################
postpred.SRMR.p <- mean(postpred.values > realized.values)


#############################################################################################################################
# Write out the p-values
#############################################################################################################################
file.name <- "SRMR.ppp.value.out"
write.table(postpred.SRMR.p, file=paste(PPMC.folder, file.name, sep=""), row.names=FALSE, col.names=FALSE)



