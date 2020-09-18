#################################################################################################################################################
# R code to compute the PPMC analyses for the IRT model
#################################################################################################################################################



################################################################################################# 
# Define a function to calculate the Compensatory normal-ogive MIRT model probability 
# for a matrix (nxJ)

# for a given matrix of thetas, matrix of a's, a vector of d's, and a vector of c's
################################################################################################

comp.3P.MIRT.normal.ogive.p.matrix <- function(theta.matrix, a.matrix, d.vector, c.vector){

  
  # theta.matrix is the matrix of thetas (for the persons, dimensions)
  # a is the matrix of discrimination parameters (for the items, dimensions)
  # d is the vector of location parameters (for the items)
  # c is the vector of lower asymptotes (for the items)
  n.local <- nrow(theta.matrix)
  J.local <- nrow(a.matrix)
  
  ones.vector <- matrix(1, nrow=n.local, ncol=1)
  
	if(ncol(theta.matrix) != ncol(a.matrix)) print("Warning: Theta and Discrimination matrices not compatible") 	
	
  phi <- pnorm(theta.matrix %*% t(a.matrix) + ones.vector %*% t(d.vector))
  
  p <- (ones.vector %*% t(c.vector)) + (1-(ones.vector %*% t(c.vector)))*phi
  
  p

}
################################################################################################




################################################################################################# 
# Define a function to calculate the Compensatory normal-ogive MIRT model probability 
# for a given vector of thetas, vector of a's, a d, and a c
# This function will be called later to evaluate the probability that any given person
# gets any given item correct
################################################################################################

comp.3P.MIRT.normal.ogive.p <- function(theta.i, a.j, d.j, c.j){

  # theta.i is the vector of thetas (for the person)
  # a.j is the vector of discrimination parameters (for the item)
	# d.j is the location parameter (for the item)
  # c.j is the lower asymptote (for the item)

	if(length(theta.i) != length(a.j)) print("Warning: Theta and Discrimination vectors of different size") 	
	
	p <- c.j+((1-c.j)*pnorm(sum(a.j*theta.i)+d.j, mean=0, sd=1, lower.tail=TRUE, log.p=FALSE))
	p

}
################################################################################################



#################################################################################################################################################

#      		GENERATE POSTERIOR PREDICTIVE DATASETS
#         & MODEL-IMPLIED EXPECTED VALUES MATRICES

#################################################################################################################################################

############################################################################################################################
# Setup a subfolder for posterior predictive datasets
############################################################################################################################
PP.datasets.folder <- paste(PPMC.folder, "PP Datasets\\", sep="")
dir.create(PP.datasets.folder)

  
############################################################################################################################
# Loop over iterations
############################################################################################################################
for(which.iter in 1:n.iters.PPMC){

  print(paste("Generating Posterior Predicted Dataset for iteration ", which.iter, sep=""))
  
  #################################################################################################################################################
  # Extract the values of the needed model parameters for this iteration
  #################################################################################################################################################
  a <- matrix(a.iterations.for.PPMC[which.iter, , ], nrow=J)
  d <- d.iterations.for.PPMC[which.iter, ]
  c <- c.iterations.for.PPMC[which.iter, ]
  theta <- matrix(theta.iterations.for.PPMC[which.iter, , ], nrow=n)

  
  #################################################################################################################################################
  # Calculate the model-implied probabilities (of a value of 1) for each data point
  #################################################################################################################################################
  probability.matrix <- comp.3P.MIRT.normal.ogive.p.matrix(
    theta.matrix=theta, 
    a.matrix=a, 
    d.vector=d, 
    c.vector=c
  )

  
  #################################################################################################################################################
  # Define the expected values in a matrix
  #################################################################################################################################################
  expected.values <- probability.matrix
  
  #################################################################################################################################################
  # Simulate uinform(0,1) values in a matrix
  #################################################################################################################################################
  u.matrix <- matrix(runif(n*J), nrow=n, ncol=J)
  
  
  #################################################################################################################################################
  # Construct the posterior predicted dataset
  #################################################################################################################################################
  x.postpred <- matrix(as.numeric(probability.matrix > u.matrix), ncol=J, byrow=FALSE)


  
  #################################################################################################################################################
  # Write out the 
  #   posterior predicted data
  #   model-implied expected values
  #################################################################################################################################################
  write.table(x.postpred, file=paste(PP.datasets.folder, "x.postpred.", which.iter, ".dat", sep=""), row.names=FALSE, col.names=FALSE)	
	write.table(expected.values, file=paste(PP.datasets.folder, "expected.values.", which.iter, ".dat", sep=""), row.names=FALSE, col.names=FALSE)	

  
} # closes loop over iterations
  
  




#################################################################################################################################################

#      		ANALYSIS OF SGDDM & SMBC

#################################################################################################################################################



################################################################################################# 
# Define a function to compute SGDDM 
# It requires
#  data.matrix
#	expected.value.matrix
################################################################################################

calculate.SGDDM <- function(data.matrix, expected.value.matrix){
	
	J.local = ncol(data.matrix)

	SMBC.matrix <- calculate.SMBC.matrix(data.matrix, expected.value.matrix)
	
	SGDDM = sum(abs((lower.tri(SMBC.matrix, diag=FALSE))*SMBC.matrix))/((J.local*(J.local-1))/2)

	SGDDM

} # closes calculate.SGDDM
################################################################################################


################################################################################################# 
# Define a function to compute a matrix of Standaridized Model-Based Covariances  
# It requires
#  data.matrix
#	expected.value.matrix
################################################################################################

calculate.SMBC.matrix <- function(data.matrix, expected.value.matrix){
	
	N.local <- nrow(data.matrix)

	MBC.matrix <- (t(data.matrix-expected.value.matrix) %*% (data.matrix-expected.value.matrix))/N.local

	MBStddevs.matrix <- diag(sqrt(diag(MBC.matrix)))

	#SMBC.matrix <- solve(MBStddevs.matrix) %*% MBC.matrix %*% solve(MBStddevs.matrix)


	J.local <- ncol(data.matrix)

	SMBC.matrix <- matrix(NA, nrow=J.local, ncol=J.local)

	for(j in 1:J.local){
		for(jj in 1:J.local){
			SMBC.matrix[j,jj] <- MBC.matrix[j,jj]/(MBStddevs.matrix[j,j]*MBStddevs.matrix[jj,jj])
		}
	}

	SMBC.matrix 

} # closes calculate.MBC.matrix
################################################################################################



########################################################################################################################
# Set up arrays to store realized and posterior predicted
#  SMBC (as a matrix)
#  SGDDM (scalar)
########################################################################################################################
realized.SMBC.array <- array(NA, c(n.iters.PPMC, J, J))
postpred.SMBC.array <- array(NA, c(n.iters.PPMC, J, J))
realized.SGDDM.vector <- array(NA, c(n.iters.PPMC))
postpred.SGDDM.vector <- array(NA, c(n.iters.PPMC))


############################################################################################################################
# Loop over iterations
############################################################################################################################
for(which.iter in 1:n.iters.PPMC){

  print(paste("Conducting PPMC for SGDDM and SMBC for iteration ", which.iter, sep=""))
  
  ############################################################################################################################
  # Read in the posterior predicted dataset
  ############################################################################################################################
  x.postpred <- as.matrix(read.table(file=paste(PP.datasets.folder, "x.postpred.", which.iter, ".dat", sep="")))
  expected.values <- as.matrix(read.table(file=paste(PP.datasets.folder, "expected.values.", which.iter, ".dat", sep=""))	)

  
  
  ########################################################################################################################
	# 						ANALYSIS OF SMBC
	########################################################################################################################

	########################################################################################################################
	# Calculate and store realized and posterior predictive values of SMBC
	########################################################################################################################

	realized.SMBC.matrix <- calculate.SMBC.matrix(
					data.matrix = x, 
					expected.value.matrix = expected.values
	)

	postpred.SMBC.matrix <- calculate.SMBC.matrix(
					data.matrix = x.postpred, 
					expected.value.matrix = expected.values
	)
  
  realized.SMBC.array[which.iter,,] <- realized.SMBC.matrix
  postpred.SMBC.array[which.iter,,] <- postpred.SMBC.matrix
  
  
  
  ########################################################################################################################
	# 						ANALYSIS OF SGDDM 
	########################################################################################################################

	########################################################################################################################
	# Calculate and store the realized and posterior predictive values of SGDDM 
	########################################################################################################################

	realized.SGDDM <- calculate.SGDDM(
					data.matrix = x, 
					expected.value.matrix = expected.values
	)

	postpred.SGDDM <- calculate.SGDDM(
					data.matrix = x.postpred, 
					expected.value.matrix = expected.values
	)
  
  realized.SGDDM.vector[which.iter] <- realized.SGDDM
  postpred.SGDDM.vector[which.iter] <- postpred.SGDDM

  
} # closes loop over iterations

  ########################################################################################################################
	# Calculate the p values for 
  #   SMBC
  #   SGDDM
	########################################################################################################################
  p.value.SMBC.matrix <- matrix(NA, nrow=J, ncol=J)
  temp <- (postpred.SMBC.array > realized.SMBC.array)
  dim(temp)
  for(j in 1:J){
    for(jj in 1:J){
      if(j != jj){
        p.value.SMBC.matrix[j,jj] <- mean(temp[ ,j,jj])
        p.value.SMBC.matrix[jj,j] <- mean(temp[ ,j,jj])        
      }
    }
  }
  
  p.value.SGDDM = mean(postpred.SGDDM.vector > realized.SGDDM.vector)

 
#############################################################################################
# Write out 
#   Summary of realized values
#   Plot
#############################################################################################


realized.values <- realized.SGDDM.vector
postpred.values <- postpred.SGDDM.vector



#############################################################################################################################
# Summarize the realized values of SGDDM
#############################################################################################################################
(SGDDM.summary <- c(summary(realized.values), HPDinterval(as.mcmc(realized.values))))

SGDDM.summary <- matrix(c(summary(realized.values), HPDinterval(as.mcmc(realized.values))), nrow=1)
colnames(SGDDM.summary)= c(
  "minimum",
  "25_percentile",
  "median",
  "mean", 
  "75_percentile",
  "maximum",
  "2.5_percentile",
  "97.5_percentile"
)  

#############################################################################################################################
# Write out the summary
#############################################################################################################################
file.name <- "Realized.SGDDM.summary.out"
write.table(SGDDM.summary, file=paste(PPMC.folder, file.name, sep=""), row.names=FALSE, col.names=TRUE, quote=FALSE)


#################################################################################################################################################
# PLOT RESULTS
#################################################################################################################################################

###########################################################################
#  Define features of the plots 
###########################################################################
plot.width=3
plot.height=3
plot.margins=c(1.5,1.5,.5,.5)+.1
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
file.main.name <- paste("Scatterplot of SGDDM", sep="")




###########################################################################
# Open a window, for the plot, which will then be saved as file types
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
  mtext(side = 1, line = 1.5, text=bquote("Realized SGDDM"))
  mtext(side = 2, line = 1.5, text=bquote("Posterior Predicted SGDDM"))
} # closes switch for if writing axis labels


#############################################################################################################################
# Copy to a .pdf
#############################################################################################################################
file.name <- paste(output.folder, file.main.name, ".pdf", sep="")
dev.copy2pdf(file=file.name)


dev.off()  # closes the device opened in R








#################################################################################################################################################
# PLOT RESULTS FOR THE SMBCs
#################################################################################################################################################
realized.values <- realized.SMBC.array
postpred.values <- postpred.SMBC.array

for(j in 1:J){
  realized.values[ ,j, j] <- NA
  postpred.values[ ,j, j] <- NA
}


###########################################################################
# Define features of the plot
###########################################################################
plot.width=6
plot.height=6
line.width=1

plot.margins=c(2, 2, .1, .1)

point.size=.6
point.type=19

###########################################################################
# Define the main name of the output files (plots of various file types)
###########################################################################
file.main.name <- paste("Plot of SMBCs (below)", sep="")



###########################################################################
# Open a window, for the plot, which will then be saved as file types
###########################################################################
windows(width=plot.width, height=plot.height)
par(family="serif")
par(mar=plot.margins)

mat=matrix(seq(1:(J*J)), ncol=J)
layout(mat=mat, 
       widths=rep(1, ncol(mat)),
       heights=rep(1, nrow(mat)), 
       respect=FALSE
)


for(j in 1:J){
  for(jj in 1:J){
    
    
    
    #############################################################################################################################
    # Plot nothing above the main diagonal
    #############################################################################################################################
    if(j>jj){
      
      #############################################################################################################################
      # Plot realized and posterior predicted values
      #############################################################################################################################
      x.values.to.plot <- realized.values[ ,j, jj]
      y.values.to.plot <- postpred.values[ ,j, jj]
      plot(
        x=x.values.to.plot,
        y=y.values.to.plot,
        ann=FALSE,
        axes=FALSE,
        type="n"
      )
      
    }    
    
    
    #############################################################################################################################
    # Plot the values above the main diagonal
    #############################################################################################################################
    if(j<jj){
      
      #############################################################################################################################
      # Plot realized and posterior predicted values
      #############################################################################################################################
      x.values.to.plot <- realized.values[ ,j, jj]
      y.values.to.plot <- postpred.values[ ,j, jj]
      #vertical.line.to.plot <- cor(X)[j,jj]
      plot(
        x=x.values.to.plot,
        y=y.values.to.plot,
        ann=FALSE,
        axes=FALSE,
        xlim=c(
          min(floor(100*min(realized.values, na.rm=TRUE))/100, floor(100*min(postpred.values, na.rm=TRUE))/100),
          max(ceiling(100*max(realized.values, na.rm=TRUE))/100, ceiling(100*max(postpred.values, na.rm=TRUE))/100)
        ),
        ylim=c(
          min(floor(100*min(realized.values, na.rm=TRUE))/100, floor(100*min(postpred.values, na.rm=TRUE))/100),
          max(ceiling(100*max(realized.values, na.rm=TRUE))/100, ceiling(100*max(postpred.values, na.rm=TRUE))/100)
        ),
        cex=point.size,
        pch=point.type,
      )
      box()
      if(jj==J) {
        axis(side=1, font=6, tck=-.03, labels=NA)
        axis(side = 1, lwd = 0, line = -.8, font=6, las=1)
      }
      if(j==1){
        axis(side=2, font=6, tck=-.03, labels=NA, las=1)
        axis(side = 2, lwd = 0, line = -.8, font=6, las = 1)
      }

      #############################################################################################################################
      # add unit line 
      #############################################################################################################################
      abline(a=0, b=1)
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
      text(x=0, y=0, var.labels[j], font=6, cex=1.5)
    }
    
  } # closes first loop over plotting 
} # closes second loop over plotting 


#############################################################################################################################
# Copy to a .pdf
#############################################################################################################################
file.name <- paste(output.folder, file.main.name, ".pdf", sep="")
dev.copy2pdf(file=file.name)


dev.off()  # closes the device opened in R




#################################################################################################################################################

#          ANALYSIS OF RAW SCORE DISTRIBUTION

#################################################################################################################################################

values <- seq(0,J,1)
raw.score.frequency.postpred.table <- NULL


############################################################################################################################
# Loop over iterations
############################################################################################################################
for(which.iter in 1:n.iters.PPMC){

  
  print(paste("Conducting PPMC for Raw Score Distribution for iteration ", which.iter, sep=""))

  
  ############################################################################################################################
  # Read in the posterior predicted dataset
  ############################################################################################################################
  x.postpred <- read.table(file=paste(PP.datasets.folder, "x.postpred.", which.iter, ".dat", sep=""))


  ############################################################################################################################
  # Compute the number correct
  ############################################################################################################################
  raw.score.postpred <- rowSums(x.postpred)

  ############################################################################################################################
  # Compute frequencies of possible scores
  ############################################################################################################################
  hist <- hist(raw.score.postpred, plot=FALSE, breaks=seq(min(values)-.5, max(values)+.5, 1))
  
  
  ############################################################################################################################
  # Make table with iteration, value, and frequency
  ############################################################################################################################
  temp <- cbind(
    rep(which.iter, length(values)), 
    values,
    hist$counts
  )
  
  raw.score.frequency.postpred.table <- rbind(raw.score.frequency.postpred.table, temp)
  
} # closes loop over iterations
  
############################################################################################################################
# Define the names of the columns of the table
############################################################################################################################
colnames(raw.score.frequency.postpred.table) <- c(
  "iteration",
  "raw.score",
  "frequency"
)  

############################################################################################################################
# Convert to a data.frame
############################################################################################################################
raw.score.frequency.postpred.table <- as.data.frame(raw.score.frequency.postpred.table)

############################################################################################################################
# Compute number correct 
# and frequencies
#   for realized data
############################################################################################################################

raw.score.realized <- rowSums(x)
frequencies.realized <- hist(raw.score.realized, plot=FALSE, breaks=seq(min(values)-.5, max(values)+.5, 1))$counts



############################################################################################################################
# Make boxplot
############################################################################################################################



###########################################################################
#  Define features of the plots 
###########################################################################
plot.width=3.5
plot.height=3
plot.margins=c(2,2.5,0,0)+.1
if(write.axis.labels==TRUE){
  plot.margins=plot.margins+c(1,1,0,0)
} # closes switch for if writing axis labels

point.size=1.5
point.type=19
xlimit<-range(realized.values,postpred.values)
ylimit<-xlimit




###########################################################################
#  Define main name of output file
###########################################################################
file.main.name <- paste("Plot of frequency of total score", sep="")


###########################################################################
# Open a window, for creating the plot, which will then be saved as file types
###########################################################################
windows(width=plot.width, height=plot.height)
par(family="serif")
par(mar=plot.margins)

  ############################################################################################################################
  # Grab the elements of the box plot
  ############################################################################################################################
  z <- boxplot(
    formula=frequency~raw.score,
    data=raw.score.frequency.postpred.table,
    plot=FALSE,
    names=as.character(values),
    xlab="Raw Score",
    ylab="Number of Examinees"
   )

  ############################################################################################################################
  # Modify so that it plots the percentiles: 
  # 2.5, 25, 50, 75, 97.5
  ############################################################################################################################

  for(which.value in 1:ncol(z$stats)){
    z$stats[,which.value] <- quantile(
      subset(raw.score.frequency.postpred.table, raw.score==values[which.value])[["frequency"]],
      probs=c(.025, .25, .5, .75, .975)
    )
    
  }

  ############################################################################################################################
  # Now plot it
  ############################################################################################################################
  
  bxp(
      z,
      outline=FALSE,
      names=as.character(values),    
      axes=FALSE,
      
    )

############################################################################################################################
# Add line for realized value
############################################################################################################################
  
lines(x=values+1, y=frequencies.realized, type="b", pch=19, cex=1.5)


  box()
  axis(side=1, font=6, tck=-.05, labels=NA)
  axis(side = 1, lwd = 0, line = -.2, font=6, at=values+1, labels=as.character(values), cex.axis=1)
  mtext(side = 1, font=6, line = 2, text="Raw Score")
  axis(side=2, font=6, tck=-.05, labels=NA, las=1)
  axis(side = 2, lwd = 0, line = -.2, font=6, las = 1, cex.axis=1)
  mtext(side = 2, font=6, line = 2.5, text="Number of Examinees")


#############################################################################################################################
# Copy to a .pdf
#############################################################################################################################
file.name <- paste(output.folder, file.main.name, ".pdf", sep="")
dev.copy2pdf(file=file.name)


dev.off()  # closes the device opened in R





#################################################################################################################################################

#          ANALYSIS OF ITEM FIT 

#################################################################################################################################################

values <- seq(0,J,1)
raw.score.frequency.postpred.table <- NULL


############################################################################################################################
# Loop over iterations
############################################################################################################################

for(which.iter in 1:n.iters.PPMC){

  
  print(paste("Conducting PPMC for Item Fit for iteration ", which.iter, sep=""))

  
  ############################################################################################################################
  # Read in the posterior predicted dataset
  ############################################################################################################################
  x.postpred <- read.table(file=paste(PP.datasets.folder, "x.postpred.", which.iter, ".dat", sep=""))


  ############################################################################################################################
  # Compute the number correct
  ############################################################################################################################
  raw.score.postpred <- rowSums(x.postpred)

  ############################################################################################################################
  # For each value of number correct
  # Compute the proportion correct for each item
  ############################################################################################################################
  temp <- NULL
  which.value=1
  for(which.value in 1:length(values)){
    current.value <- values[which.value]
    x.postpred.current.value <- subset(x.postpred, raw.score.postpred==current.value)
    temp <- rbind(temp, c(which.value-1, colMeans(x.postpred.current.value)))
    
  }
 

  ############################################################################################################################
  # Make table with iteration, value, and proportions
  ############################################################################################################################
  temp <- cbind(rep(which.iter, length(values)), temp)
  
  ############################################################################################################################
  # Add it to the previous values
  ############################################################################################################################
  raw.score.frequency.postpred.table <- rbind(raw.score.frequency.postpred.table, temp)
  
} # closes loop over iterations
  

############################################################################################################################
# Define the names of the columns of the table
############################################################################################################################
colnames(raw.score.frequency.postpred.table) <- c(
  "iteration",
  "raw.score",
  paste("proportion_correct_item_", seq(1,J), sep="")  
)  


############################################################################################################################
# Convert to a data.frame
############################################################################################################################
raw.score.frequency.postpred.table <- as.data.frame(raw.score.frequency.postpred.table)

############################################################################################################################
# Compute number correct 
# and proportions correct
#   for realized data
############################################################################################################################

raw.score.realized <- rowSums(x)
raw.score.frequency.realized.table <- NULL
which.value=1
for(which.value in 1:length(values)){
  current.value <- values[which.value]
  x.realized.current.value <- subset(x, raw.score.realized==current.value)
  raw.score.frequency.realized.table <- rbind(raw.score.frequency.realized.table, c(current.value, colMeans(x.realized.current.value)))
  
}
colnames(raw.score.frequency.realized.table) <- c(
  "raw.score", 
  paste("proportion_correct_item_", seq(1,J), sep="")  
  )

raw.score.frequency.realized.table <- as.data.frame(raw.score.frequency.realized.table)





############################################################################################################################
# Make boxplot
############################################################################################################################

###########################################################################
#  Define features of the plots 
###########################################################################
plot.width=6.5
plot.height=3.5
plot.margins=c(3.5,5.5,2.0,0)+.1


point.size=1.5
point.type=19
caxis=1.5


###########################################################################
#  Define main name of output file
###########################################################################
file.main.name <- paste("Plot of proportion corrects by total score", sep="")


###########################################################################
# Open a window, for creating the plot, which will then be saved as file types
###########################################################################
windows(width=plot.width, height=plot.height)
par(family="serif")
par(mar=plot.margins)
  
  


mat=matrix(seq(1:(J+1)), ncol=3, byrow=TRUE)
layout(mat=mat, 
       respect=FALSE
)

j=1
for(j in 1:J){
  ############################################################################################################################
  # Grab the elements of the box plot
  ############################################################################################################################
  z <- boxplot(
    formula=as.formula(paste("proportion_correct_item_", j, "~raw.score", sep="")),
    data=raw.score.frequency.postpred.table,
    plot=FALSE,
    names=as.character(values),
    axes=FALSE
  )
  

  ############################################################################################################################
  # Modify so that it plots the percentiles: 
  # 2.5, 25, 50, 75, 97.5
  ############################################################################################################################
  which.value=1
  for(which.value in 1:ncol(z$stats)){
    z$stats[,which.value] <- quantile(
      # what=subset(raw.score.frequency.postpred.table, raw.score==values[which.value])
      subset(raw.score.frequency.postpred.table, raw.score==values[which.value])[[(paste("proportion_correct_item_", j , sep=""))]],
      probs=c(.025, .25, .5, .75, .975),
      na.rm=TRUE
    )
    
  }

  ############################################################################################################################
  # Now plot it
  ############################################################################################################################
  
  bxp(
      z,
      outline=FALSE,
      names=as.character(values),    
      #font.axis=6,
      axes=FALSE,      
    )
  
  ############################################################################################################################
  # Add line for realized value
  ############################################################################################################################
  lines(x=values+1, y=raw.score.frequency.realized.table[[paste("proportion_correct_item_", j, sep="")]], type="b", pch=point.type, cex=point.size)

  box()
  axis(side=1, font=6, tck=-.06, labels=NA)
  
  axis(side = 1, lwd = 0, line = -.2, font=6, at=values+1, labels=as.character(values), cex.axis=caxis)
  mtext(side = 1, font=6, line = 2, text="Raw Score")
  axis(side=2, font=6, tck=-.06, labels=NA, las=1)
  axis(side = 2, lwd = 0, line = -.2, font=6, las = 1, cex.axis=caxis)
  mtext(side = 2, font=6, line = 3, text="Prop. Correct")
  mtext(side = 3, font=6, line = .5, text=paste("Item ", j, sep=""))
 
  
  
} # closes loop over j
 
#############################################################################################################################
# Copy to a .pdf
#############################################################################################################################
file.name <- paste(output.folder, file.main.name, ".pdf", sep="")
dev.copy2pdf(file=file.name)


dev.off()  # closes the device opened in R



