
################################################################
## Function to compute CIs for log-rates and e0
################################################################

sim_to_ci <- function(ETAsim,E0sim,alpha){
  ## dimensions
  m <- dim(ETAsim)[1]
  nc <- dim(ETAsim)[2]
  ## alpha level
  alpha.lev <- alpha/100
  ## median, lower and upper ETA
  ETA.sim.upp <- ETA.sim.low <- ETA.sim.med <- matrix(NA,m,nc)
  i <- 1
  for (i in 1:m){
    ETA.sim.med[i,] <- apply(ETAsim[i,,], 1, median, na.rm = T)
    ETA.sim.upp[i,] <- apply(ETAsim[i,,], 1, quantile, prob=1- (1-alpha.lev)/2, na.rm = T)
    ETA.sim.low[i,] <- apply(ETAsim[i,,], 1, quantile, prob=(1-alpha.lev)/2, na.rm = T)
  }
  ## focus on e0 only
  E0sim.age0 <- E0sim[1,,]
  e0.sim.med <- apply(E0sim.age0, 1, median, na.rm = T)
  e0.sim.upp <- apply(E0sim.age0, 1, quantile, prob=1- (1-alpha.lev)/2, na.rm = T)
  e0.sim.low <- apply(E0sim.age0, 1, quantile, prob=(1-alpha.lev)/2, na.rm = T)
  ## output
  out <- list(ETA.sim.med=ETA.sim.med,ETA.sim.upp=ETA.sim.upp,ETA.sim.low=ETA.sim.low,
              e0.sim.med=e0.sim.med,e0.sim.upp=e0.sim.upp,e0.sim.low=e0.sim.low)
  return(out)
}