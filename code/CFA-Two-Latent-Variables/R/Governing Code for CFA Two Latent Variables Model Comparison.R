########################################################################################################################
# Governing R code to run model comparison for CFA example with two latent variables

# data file should be in the Data folder
# coda files for 3 chains should be in the coda folder
########################################################################################################################

rm(list=ls())



########################################################################################################################
# Define the 
#   main analysis folder
#   R folder
#   Data folder
#   Coda folder
########################################################################################################################
main.folder <- "C:\\BPM\\CFA Two Latent Variables\\"
R.folder <- paste(main.folder, "R\\", sep="")
data.folder <- paste(main.folder, "Data\\", sep="")
coda.folder <- paste(main.folder, "coda files\\", sep="")



########################################################################################################################
# Read in the observed data
########################################################################################################################
file.name <- "IIS.dat"
x <- as.matrix(read.table(paste(data.folder, file.name, sep=""), header=TRUE))




############################################################################################################################
# Extract features
# n is the number of examinees
# J is the number of observables
############################################################################################################################
n <- nrow(x)
J <- ncol(x)





########################################################################################################################
# Conduct Bayes Factor analyses
#   Create a subfolder if one isn't there 
########################################################################################################################
Bayes.factor.folder <- paste(main.folder, "Bayes Factor\\", sep="")
dir.create(Bayes.factor.folder)
setwd(Bayes.factor.folder)


  ###########################################################################
  #  Read in the draws
  #  Using the read.coda function as the index name is different
  ###########################################################################
  setwd(coda.folder)
  library(coda)
  chain.1 <- read.coda("codaforCPO1.txt", "codaforCPOIndex.txt")
  chain.2 <- read.coda("codaforCPO2.txt", "codaforCPOIndex.txt")
  chain.3 <- read.coda("codaforCPO3.txt", "codaforCPOIndex.txt")
  draws.from.bugs <- mcmc.list(chain.1, chain.2, chain.3)
  
  #########################################################################################################################
  # Convert the draws to a matrix
  #########################################################################################################################
  draws.from.bugs.as.matrix <- as.matrix(draws.from.bugs)
  
  
  
  ###########################################################################
  #  Declare burnin and thin if needed
  ###########################################################################
  n.burnin=500
  
  
  #############################################################################################################################
  # Select the iterations to summarize
  #############################################################################################################################
  draws.to.analyze <- window(draws.from.bugs,
    start=n.burnin+1)
  
  #########################################################################################################################
  # Convert the draws to a matrix
  #########################################################################################################################
  draws.to.analyze.as.matrix <- as.matrix(draws.to.analyze)
  dim(draws.to.analyze.as.matrix)
  
  
  ###########################################################################
  #  Summarize the Posterior
  ###########################################################################
  file.name <- "Summarize Posterior.R"
  source(paste(R.folder, file.name, sep=""))
  
  
  ###########################################################################
  #  Compute CPO and Psuedomarginal Likelihood
  ###########################################################################
  setwd(Bayes.factor.folder)
  file.name <- "Compute CPO and Psuedomarginal Likelihood.R"
  source(paste(R.folder, file.name, sep=""))


