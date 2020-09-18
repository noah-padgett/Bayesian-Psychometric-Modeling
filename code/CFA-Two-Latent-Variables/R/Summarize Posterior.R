#################################################################################################################################################
# Source code for running normal model with BUGS 
# This code summarizes the draws from BUGS
#################################################################################################################################################


#########################################################################################################################
# Call the 'coda' library
#########################################################################################################################
library(coda)



#############################################################################################################################
# Combine chains for summaries
#############################################################################################################################

coda.options(combine.stats=TRUE, combine.plots=TRUE)


#########################################################################################################################
# Extract the summary statistics
#   Usual
#   Percentiles
#   HPD
#########################################################################################################################
summary.stats <- summary(draws.to.analyze)
summary.stats$statistics
summary.stats$quantiles

probability.for.HPD=.95

summary.statistics <- cbind(
  summary.stats$statistics, 
  summary.stats$quantiles, 
  matrix(HPDinterval(draws.to.analyze, prob=probability.for.HPD)[[1]], ncol=2)
)

colnames(summary.statistics) <- c(
  colnames(summary.stats$statistics),
  colnames(summary.stats$quantiles),
  c("95% HPD lower", "95% HPD Upper")
)

