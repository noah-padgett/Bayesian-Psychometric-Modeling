# Useful R functions for ploting and inspecting results

# for Stan, extract results and compute Gelman-Rubin-Brooks convergence plot
plot_rhat <- function(fit, par){
  dat <- as.array(fit)
  dim.dat <- dim(dat)
  rhat_est <- data.frame(matrix(nrow=dim.dat[1], ncol=dim.dat[3]+1))
  colnames(rhat_est) <- c('iter', dimnames(dat)[['parameters']])
  rhat_est$iter <- 1:dim.dat[1]

  pdat <- dat[,, parameters=par]

  for(i in 1:dim.dat[1]){
    rhat_est[i, par] <- rstan::Rhat(pdat[1:i,])
  }

  rhat_est$PARA <- rhat_est[, par]

  a<-ggplot(rhat_est, aes(x=iter, y=PARA)) +
    geom_line()+
    geom_hline(yintercept = 1, linetype="dashed", color="grey")+
    lims(y=c(min(rhat_est[,par])-0.05, max(rhat_est[,par])+0.05)) +
    labs(y="Rhat",
         x="Post-warm-up iteration",
         title=paste0("Iteration Evolution of Gelman-Rubin-Brooks Convergence Criterion: ", par),
         subtitle = "Should approach 1 as iterations increase") +
    theme_bw() +
    theme(panel.grid = element_blank())
  return(a)
}
