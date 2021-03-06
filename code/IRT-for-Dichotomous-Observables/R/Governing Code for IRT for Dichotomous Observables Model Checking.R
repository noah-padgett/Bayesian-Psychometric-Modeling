########################################################################################################################
# Governing R code to run model checking IRT examples for Dichtomous Data

# data file should be in the Data folder
# coda files for 3 chains should be in the coda folder
########################################################################################################################

rm(list=ls())


########################################################################################################################
# Define the
#   main analysis folder
#   R folder
#   Data folder
#   coda folder
########################################################################################################################
main.folder <- paste0(getwd(),"/code/IRT-for-Dichotomous-Observables/")
R.folder <- paste(main.folder, "R/", sep="")
data.folder <- paste(main.folder, "Data/", sep="")
coda.folder <- paste(main.folder, "coda files/", sep="")




########################################################################################################################
# Read in the observed data
########################################################################################################################
file.name <- "LSAT.dat"
x <- as.matrix(read.table(paste(data.folder, file.name, sep=""), header=TRUE))



############################################################################################################################
# Extract features
# n is the number of examinees
# J is the number of observables
############################################################################################################################
n <- nrow(x)
J <- ncol(x)



########################################################################################################################
# PPMC Analyses
########################################################################################################################

  ###########################################################################
  # Read in the draws from the coda folder
  ###########################################################################
  bugs.sim <-  c("coda1.txt", "coda2.txt", "coda3.txt")
  setwd(coda.folder)
  file.name <- "Read in Draws.R"
  source(paste(R.folder, file.name, sep=""))

  ###########################################################################
  # Create PPMC folder if it isn't there
  ###########################################################################
  PPMC.folder <- paste(main.folder, "PPMC Analyses\\", sep="")
  dir.create(PPMC.folder)


  ###########################################################################
  #  Declare burnin and thin if needed
  ###########################################################################
  n.burnin=6000


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

  #########################################################################################################################
  # Remove the original draws from BUGS, so as to avoid any errors
  #########################################################################################################################
  rm(draws.from.bugs, draws.from.bugs.as.matrix)



  ############################################################################################################################
  # Define the total number of iterations
  ############################################################################################################################
  n.iters.total <- nrow(draws.to.analyze.as.matrix)


  ########################################################################################################################
  # Structure the draws
  ########################################################################################################################
  file.name <- "Structure the Draws.R"
  source(paste(R.folder, file.name, sep=""))


  ########################################################################################################################
  # Thin the draws
  ########################################################################################################################
  n.thin.for.PPMC=4
  iterations.to.use.PPMC <- seq(from=1, to=n.iters.total, by=n.thin.for.PPMC)
  n.iters.PPMC <- length(iterations.to.use.PPMC)

  file.name <- "Thin the Draws.R"
  source(paste(R.folder, file.name, sep=""))



  ########################################################################################################################
  # Conduct the PPMC analyses
  #   Create a subfolder if one isn't there
  ########################################################################################################################
  setwd(PPMC.folder)
  output.folder <- PPMC.folder
  write.axis.labels=TRUE
  file.name <- "PPMC Analyses.R"
  source(paste(R.folder, file.name, sep=""))





########################################################################################################################
# Conduct BILOG analyses
########################################################################################################################

########################################################################################################################
# Define the
#   BILOG folder
########################################################################################################################
BILOG.folder <- paste(LSAT.folder, "BILOG\\", sep="")
dir.create(BILOG.folder)

########################################################################################################################
# Generate out BILOG code
# For 3 PL
# Need to modify, run by hand
########################################################################################################################
setwd(BILOG.folder)
file.name <- "3PL BILOG via irtoys.R"
source(paste(R.folder, file.name, sep=""))
