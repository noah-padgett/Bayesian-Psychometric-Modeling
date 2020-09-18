#################################################################################################################################################
# This code thins the draws 
#################################################################################################################################################


############################################################################################################################
# Define the iterations to use in PPMC
############################################################################################################################
a.iterations.for.PPMC <- array(a.iterations[iterations.to.use.PPMC,,], c(n.iters.PPMC,J,M))
d.iterations.for.PPMC <- array(d.iterations[iterations.to.use.PPMC,], c(n.iters.PPMC,J))
c.iterations.for.PPMC <- array(c.iterations[iterations.to.use.PPMC,], c(n.iters.PPMC,J))
theta.iterations.for.PPMC <- array(theta.iterations[iterations.to.use.PPMC,,], c(n.iters.PPMC,n,M))

