#################################################################################################################################################
# This code structures the draws according to the model
#################################################################################################################################################


############################################################################################################################
# Define the amount of thinning and the iterations to keep
############################################################################################################################
#n.thin.for.PPMC=4
#iterations.to.use.PPMC <- seq(from=1, to=n.iters.total, by=n.thin.for.PPMC)
#n.iters.PPMC <- length(iterations.to.use.PPMC)

############################################################################################################################
# Define the iterations to use in PPMC
############################################################################################################################
dim(gamma.iterations)
dim(pi.iterations)
dim(theta.iterations)

gamma.iterations.for.PPMC <- gamma.iterations[iterations.to.use.PPMC,]
dim(gamma.iterations.for.PPMC)

pi.iterations.for.PPMC <- pi.iterations[iterations.to.use.PPMC,,]
dim(pi.iterations.for.PPMC)

theta.iterations.for.PPMC <- theta.iterations[iterations.to.use.PPMC,]
dim(theta.iterations.for.PPMC)  
  
  