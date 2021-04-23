# model code
jags.model.cfa <- function(){

  ###################
  # Specify the factor analysis measurement model for the observables
  ####################
  for (i in 1:n){

    # expected value for each examinee for each LRV
    mu[i,1] <- lambda[1,1]*ksi[i,1]
    mu[i,2] <- lambda[2,1]*ksi[i,1]
    mu[i,3] <- lambda[3,1]*ksi[i,1]
    mu[i,4] <- lambda[4,2]*ksi[i,2]
    mu[i,5] <- lambda[5,2]*ksi[i,2]

    xstar[i, 1:J] ~ dmnorm(mu[i, 1:J], inv.psi[1:J,1:J])    # distribution for LRVs
    # obtain prob of each response
    for(j in 1:J){
      probit(p[i, j]) <- xstar[i, j] - tau[j]
      x[i,j] ~ dbern(p[i, j])
    }

  }


  ######################################################################
  # Specify the (prior) distribution for the latent variables
  ######################################################################
  for (i in 1:n){
    ksi[i, 1:M] ~ dmnorm(kappa[], inv.phi[,])  # distribution for the latent variables
  }


  ######################################################################
  # Specify the prior distribution for the parameters that govern the latent variables
  ######################################################################
  for(m in 1:M){
    kappa[m] = 0              # Means of latent variables
  }

  rho ~ dunif(-1,1)
  phi0[1,2] = rho
  phi0[2,1] = rho
  phi0[1,1] = 1
  phi0[2,2] = 1
  inv.phi[1:2,1:2] = inverse(phi0[1:2,1:2])

  ######################################################################
  # Specify the prior distribution for the measurement model parameters
  ######################################################################
  inv.psi[1:J,1:J] ~ dwish(psi.0[ , ], 6)
  psi[1:J,1:J] <- inverse(inv.psi[ , ])

  for(j in 1:J){
    tau[j] ~ dnorm(0, .1)        # thresholds for observables
    inv.psi[j,j] = 1
  }

  lambda[1,1] ~ dnorm(1, .1)
  lambda[4,2] ~ dnorm(1, .1)
  lambda[2,1] ~ dnorm(1, .1)
  lambda[3,1] ~ dnorm(1, .1)
  lambda[5,2] ~ dnorm(1, .1)
}
# data must be in a list
dat <- read.table("data/LSAT.dat", header=T)

A <- diag(1, nrow=5, ncol=5)
#diag(A) <- NA

mydata <- list(
  n = nrow(dat),
  J = ncol(dat),
  M = 2,
  x = as.matrix(dat),
  #phi.0 = matrix(c(2,0.6, 0.6, 2), ncol=2),
  psi.0 = A

  #d = 2
)


start_values <- list(
  list("tau"=c(0,0,0,0,0),
       lambda= structure(
         .Data= c( NA,  2.00E+00, 2.00E+00, NA, NA,
                   NA,  NA,  NA,  NA, 2.00E+00),
         .Dim=c(5, 2)),
       rho = 0,
       inv.psi=structure(
         .Data=c(NA, 0, 0, 0, 0,
                 0, NA, 0, 0, 0,
                 0, 0, NA, 0, 0,
                 0, 0, 0, NA, 0,
                 0, 0, 0, 0, NA),
         .Dim=c(5,5)))
  ,
  list(tau=c(2, 2, 2, 2, 2),
       lambda= structure(
         .Data= c( NA,  5.00E-01, 5.00E-01,  NA, NA,
                   NA,  NA,  NA,  NA, 5.00E-01),
         .Dim=c(5, 2)),
       rho = 0.25,
       inv.psi=structure(
         .Data=c(NA, 0, 0, 0, 0,
                 0, NA, 0, 0, 0,
                 0, 0, NA, 0, 0,
                 0, 0, 0, NA, 0,
                 0, 0, 0, 0, NA),
         .Dim=c(5,5)))
  ,
  list(tau=c(-2, -2, -2, -2, -2),
       lambda= structure(
         .Data= c( NA,  1.00E+00, 1.00E+00,  NA, NA,
                   NA,  NA,  NA,  NA, 1.00E+00),
         .Dim=c(5, 2)),
       rho = -0.25,
       inv.psi=structure(
         .Data=c(NA, 0, 0, 0, 0,
                 0, NA, 0, 0, 0,
                 0, 0, NA, 0, 0,
                 0, 0, 0, NA, 0,
                 0, 0, 0, 0, NA),
         .Dim=c(5,5)))
)


# vector of all parameters to save
param_save <- c("tau", "lambda", "rho", "psi")

# fit model
fit <- jags(
  model.file=jags.model.cfa,
  data=mydata,
  inits=start_values,
  parameters.to.save = param_save,
  n.iter=5000,
  n.burnin = 2500,
  n.chains = 3)

print(fit)

# extract posteriors for all chains
jags.mcmc <- as.mcmc(fit)


# convert to single data.frame for density plot
a <- colnames(as.data.frame(jags.mcmc[[1]]))
plot.data <- data.frame(as.matrix(jags.mcmc, chains=T, iters = T))
colnames(plot.data) <- c("chain", "iter", a)


plot_title <- ggtitle("Posterior distributions",
                      "with medians and 80% intervals")
mcmc_areas(plot.data, regex_pars = "tau", prob = 0.8) +
  plot_title

mcmc_trace(plot.data, regex_pars = "tau")

mcmc_areas(plot.data, regex_pars = "lambda", prob = 0.8) +
  plot_title
mcmc_trace(plot.data, regex_pars = "lambda")

mcmc_areas(
  plot.data,, regex_pars = "phi",
  prob = 0.8) +
  plot_title
mcmc_trace(plot.data, regex_pars = "phi")

mcmc_areas(
  plot.data,, regex_pars = "psi",
  prob = 0.8) +
  plot_title
mcmc_trace(plot.data, regex_pars = "psi")
