
## --------------------------------------------------------- ##
##
## Functions for fitting the CCP-splines model
##
## --------------------------------------------------------- ##



## function to build up B-splines and associated bases for derivatives
BsplineGrad <- function(x, xl, xr, ndx=NULL, deg, knots=NULL){
  if(is.null(knots)){
    dx <- (xr - xl)/ndx
    knots <- seq(xl - deg * dx, xr + deg * dx, by=dx)
    knots <- round(knots, 8)
  }else{
    knots <- knots
    dx <- diff(knots)[1]
  }
  P <- outer(x, knots, MortSmooth_tpower, deg)
  n <- dim(P)[2]
  D <- diff(diag(n), diff=deg+1)/(gamma(deg+1)*dx^deg)
  B <- (-1)^(deg + 1) * P %*% t(D)
  ##
  knots1 <- knots[-c(1,length(knots))]
  P <- outer(x, knots1, MortSmooth_tpower, deg-1)
  n <- dim(P)[2]
  D <- diff(diag(n),diff=deg)/(gamma(deg)*dx^(deg-1))
  BB <- ((-1)^(deg) * P %*% t(D))/dx
  D <- diff(diag(ncol(BB) + 1))
  C <- BB %*% D
  ##
  out <- list(dx=dx, knots=knots, B=B, C=C)
}

## function for estimating two-dimensional P-splines 
## addressing infant mortality 
## for a given set of smoothing parameters lambdas
PSinfant <- function(Y, E, lambdas, WEI, ALPHAS.st=NULL, verbose=FALSE){
  require(MortalitySmooth)
  ## dimensions
  m <- nrow(Y)
  n <- ncol(Y)
  a <- 1:m
  t <- 1:n
  ## w/o age 0
  a0 <- a[-1]
  m0 <- m-1

  ## replace NA in Y and E to avoid issues later one
  E[WEI==0] <- 999
  Y[WEI==0] <- 99
  
  ## B-splines basis
  ## with infant-specialized coeff
  ## over ages w/o age 0
  a0min <- min(a0)
  a0max <- max(a0)
  nda0 <- floor(m0/5)
  dega <- 3
  BCa0 <- BsplineGrad(a0, a0min, a0max, nda0, dega)
  Ba0 <- BCa0$B
  nba0 <- ncol(Ba0)
  ## adding infant-specific basis
  Ba <- cbind(0, Ba0)
  Ba <- rbind(c(1, rep(0,nba0)), Ba)
  nba <- ncol(Ba)
  ## basis for the derivatives
  ## over ages w/o age 0
  Ca0 <- BCa0$C
  ## adding infant-specific basis
  Ca <- cbind(0, Ca0)
  Ca <- rbind(c(-1,Ba[2,1:dega+1],rep(0,nba0-dega)), Ca)
  ## over years
  tmin <- min(t)
  tmax <- max(t)
  ndt <- floor(n/5)
  degt <- 3
  BCt <- BsplineGrad(t, tmin, tmax, ndt, degt)
  Bt <- BCt$B
  nbt <- ncol(Bt)
  ## basis for the derivatives
  Ct <- BCt$C
  
  ## tensor product of B-splines for the GLAM
  Ba1 <- kronecker(matrix(1, ncol=nba, nrow=1), Ba)
  Ba2 <- kronecker(Ba, matrix(1, ncol=nba, nrow=1))
  RTBa <- Ba1 * Ba2
  Bt1 <- kronecker(matrix(1, ncol=nbt, nrow=1), Bt)
  Bt2 <- kronecker(Bt, matrix(1, ncol=nbt, nrow=1))
  RTBt <- Bt1 * Bt2
  
  ## penalty terms
  ## over ages
  Da <- diff(diag(nba), diff=2)
  Da[1,1] <- 0
  tDDa <- t(Da) %*% Da
  ## over years
  Dt <- diff(diag(nbt), diff=2)
  tDDt <- t(Dt) %*% Dt
  ## kronecker product of difference matrices
  Pa <- kronecker(diag(nbt), tDDa)
  Pt <- kronecker(tDDt, diag(nba))
  ## smoothing parameters
  lambda.a <- lambdas[1]
  lambda.t <- lambdas[2]
  P <- lambda.a * Pa + lambda.t * Pt
  
  ## starting values
  if(is.null(ALPHAS.st)){
    ETA0 <- log( (Y+0.1) / (E+0.1) )
    ## for the age-time with weights=0
    y <- c(Y)
    e <- c(E)
    wei <- c(WEI)
    aa <- rep(a, n)
    tt <- rep(t, each = m)
    fit0 <- glm(round(y) ~ aa + tt + offset(log(e)), family=poisson, weights=wei)
    ETAglm <- matrix(cbind(1,aa,tt) %*% fit0$coef, m, n)
    ETA0[WEI == 0] <- ETAglm[WEI == 0]
    ## starting coefficients
    BBa <- solve(t(Ba) %*% Ba + 1e-06*diag(nba), t(Ba))
    BBt <- solve(t(Bt) %*% Bt + 1e-06*diag(nbt), t(Bt))
    ALPHAS <- MortSmooth_BcoefB(BBa, BBt, ETA0)
  }else{
    ALPHAS <- ALPHAS.st
  }
  ## Poisson iteration within a GLAM framework
  
  ## starting deviance for convergence
  dev <- 10^6
  ## modified counts to avoid -Inf in computing deviance
  Y0 <- Y
  Y0[Y==0] <- 10^-8
  for(it in 1:20){
    ETA <- MortSmooth_BcoefB(Ba, Bt, ALPHAS) ## linear predictor
    MU <- E * exp(ETA) ## expected values
    
    ## update deviance
    dev.old <- dev
    dev <- 2 * sum(WEI * (Y * log(Y0/MU)))
    d.dev <- abs(dev - dev.old)/(0.1 + abs(dev))
    
    W <- MU ## regression weights
    WW <- WEI * W
    ZZ <- ETA + (1/MU) * (Y - MU) ## working response
    ZZ[which(WEI == 0)] <- 0
    ## scoring objects
    tBWB <- MortSmooth_BWB(RTBa, RTBt, nba, nbt, WW)
    tBWBpP <- tBWB + P
    tBWZ <- MortSmooth_BcoefB(t(Ba), t(Bt), (WW * ZZ))
    alphas <- solve(tBWBpP, c(tBWZ)) ## update coefficients
    ALPHAS <- matrix(alphas, nrow = nba) ## in matrix
    if(verbose) cat(it, d.dev, "\n")
    if(d.dev < 10^-4 & it>=4) break
  }
  ## fitted values
  ALPHAS.hat <- ALPHAS
  ETA.hat <- ETA
  MU.hat <- MU
  ## fitted derivatives over ages
  ETA.hat1a <- MortSmooth_BcoefB(Ca, Bt, ALPHAS.hat)
  ## fitted derivatives over years
  ETA.hat1t <- MortSmooth_BcoefB(Ba, Ct, ALPHAS.hat)
  ## diagnostics
  H <- solve(tBWBpP, tBWB)
  h <- diag(H)
  ed <- sum(h)
  psi <- dev / (sum(WEI) - ed)
  bic <- dev + log(sum(WEI))*ed
  
  ## return objects
  out <- list(
    ## original data
    Y=Y, E=E, lambdas=lambdas,
    ## diagnostic
    dev=dev, ed=ed, bic=bic, psi=psi,
    ## fitted values
    ALPHAS=ALPHAS.hat, ETA=ETA.hat, MU=MU.hat,
    ETA1a=ETA.hat1a, ETA1t=ETA.hat1t, 
    ## basis and associated objects
    Ba=Ba, Bt=Bt,
    Ca=Ca, Ct=Ct,
    ## convergence objects
    d.dev=d.dev, it=it
  )
  return(out)
}

## function for estimating CP-splines
## for a given set of smoothing parameters lambdas
## commonly taken from the estimated object, PSinfant()
CPSfunction <- function(Y, E, WEI, lambdas, 
                        kappas=c(10^4, 10^4),
                        deltas, S,
                        verbose=FALSE,
                        ALPHAS.st=NULL){
  ## dimensions
  m <- nrow(Y)
  n <- ncol(Y)
  a <- 1:m
  t <- 1:n
  ## w/o age 0
  a0 <- a[-1]
  m0 <- m-1
  
  ## replace NA in Y and E to avoid issues later one
  E[WEI==0] <- 999
  Y[WEI==0] <- 99
  
  ## over ages w/o age 0
  a0min <- min(a0)
  a0max <- max(a0)
  nda0 <- floor(m0/5)
  dega <- 3
  BCa0 <- BsplineGrad(a0, a0min, a0max, nda0, dega)
  Ba0 <- BCa0$B
  nba0 <- ncol(Ba0)
  ## adding infant-specific basis
  Ba <- cbind(0, Ba0)
  Ba <- rbind(c(1, rep(0,nba0)), Ba)
  nba <- ncol(Ba)
  ## basis for the derivatives
  ## over ages w/o age 0
  Ca0 <- BCa0$C
  ## adding infant-specific basis
  Ca <- cbind(0, Ca0)
  Ca <- rbind(c(-1,Ba[2,1:dega+1],rep(0,nba0-dega)), Ca)
  ## over years
  tmin <- min(t)
  tmax <- max(t)
  ndt <- floor(n/5)
  degt <- 3
  BCt <- BsplineGrad(t, tmin, tmax, ndt, degt)
  Bt <- BCt$B
  nbt <- ncol(Bt)
  ## basis for the derivatives
  Ct <- BCt$C
  
  ## tensor product of B-splines for the GLAM
  Ba1 <- kronecker(matrix(1, ncol=nba, nrow=1), Ba)
  Ba2 <- kronecker(Ba, matrix(1, ncol=nba, nrow=1))
  RTBa <- Ba1 * Ba2
  Bt1 <- kronecker(matrix(1, ncol=nbt, nrow=1), Bt)
  Bt2 <- kronecker(Bt, matrix(1, ncol=nbt, nrow=1))
  RTBt <- Bt1 * Bt2
  ## tensor product for the derivatives
  Ca1 <- kronecker(matrix(1, ncol=nba, nrow=1), Ca)
  Ca2 <- kronecker(Ca, matrix(1, ncol=nba, nrow=1))
  RTCa <- Ca1 * Ca2
  Ct1 <- kronecker(matrix(1, ncol=nbt, nrow=1), Ct)
  Ct2 <- kronecker(Ct, matrix(1, ncol=nbt, nrow=1))
  RTCt <- Ct1 * Ct2
  
  ## penalty terms for smoothing
  ## over ages
  Da <- diff(diag(nba), diff=2)
  ## no penalization over age for age 0
  Da[1,1] <- 0
  tDDa <- t(Da) %*% Da
  ## over years
  Dt <- diff(diag(nbt), diff=2)
  tDDt <- t(Dt) %*% Dt
  ## kronecker product of difference matrices
  Pa <- kronecker(diag(nbt), tDDa)
  Pt <- kronecker(tDDt, diag(nba))
  ## smoothing parameters
  lambda.a <- lambdas[1]
  lambda.t <- lambdas[2]
  P <- lambda.a * Pa + lambda.t * Pt
  
  ## penalty terms for constraints
  ## extract deltas
  delta.a.up <- deltas$delta.a.up
  delta.a.low <- deltas$delta.a.low
  delta.t.up <- deltas$delta.t.up
  delta.t.low <- deltas$delta.t.low
  ## construct G matrices
  ones12 <- matrix(1, n, 1)
  g.a.up <- kronecker(ones12, delta.a.up)
  G.a.up <- matrix(g.a.up, m, n)
  g.a.low <- kronecker(ones12, delta.a.low)
  G.a.low <- matrix(g.a.low, m, n)
  g.t.up <- kronecker(ones12, delta.t.up)
  G.t.up <- matrix(g.t.up, m, n)
  g.t.low <- kronecker(ones12, delta.t.low)
  G.t.low <- matrix(g.t.low, m, n)
  
  if(is.null(ALPHAS.st)){
    ETA0 <- log( (Y+0.1) / (E+0.1) )
    ## for the age-time with weights=0
    y <- c(Y)
    e <- c(E)
    wei <- c(WEI)
    aa <- rep(a, n)
    tt <- rep(t, each=m)
    fit0 <- glm(round(y) ~ aa + tt + offset(log(e)), family=poisson, weights=wei)
    ETAglm <- matrix(cbind(1,aa,tt) %*% fit0$coef, m, n)
    ETA0[WEI == 0] <- ETAglm[WEI == 0]
    ## starting coefficients
    BBa <- solve(t(Ba) %*% Ba + 1e-06*diag(nba), t(Ba))
    BBt <- solve(t(Bt) %*% Bt + 1e-06*diag(nbt), t(Bt))
    ALPHAS <- MortSmooth_BcoefB(BBa, BBt, ETA0)
  }else{
    ALPHAS <- ALPHAS.st
  }
  ## Poisson iteration with asymmetric penalty in a GLAM framework
  
  ## starting deviance for convergence
  dev <- 10^6
  ## modified counts to avoid -Inf in computing deviance
  Y0 <- Y
  Y0[Y==0] <- 10^-8
  
  for(it in 1:100){
    ## over ages
    CaBt.ALPHAS <- MortSmooth_BcoefB(Ca, Bt, ALPHAS)
    ## up
    v.a.up <- CaBt.ALPHAS > G.a.up
    v.a.up <- v.a.up * S
    P.a.up <- MortSmooth_BWB(RTCa, RTBt, nba, nbt, v.a.up)
    P.a.up <- kappas[1] * P.a.up
    p.a.up <- MortSmooth_BcoefB(t(Ca), t(Bt), v.a.up*G.a.up)
    p.a.up <- kappas[1] * p.a.up
    ## low
    v.a.low <- CaBt.ALPHAS < G.a.low
    v.a.low <- v.a.low * S
    P.a.low <- MortSmooth_BWB(RTCa, RTBt, nba, nbt, v.a.low)
    P.a.low <- kappas[1] * P.a.low
    p.a.low <- MortSmooth_BcoefB(t(Ca), t(Bt), v.a.low*G.a.low)
    p.a.low <- kappas[1] * p.a.low
    
    P.a <- P.a.up + P.a.low
    p.a <- p.a.up + p.a.low
    ## over years
    BaCt.ALPHAS <- MortSmooth_BcoefB(Ba, Ct, ALPHAS)
    ## up
    v.t.up <- BaCt.ALPHAS > G.t.up
    v.t.up <- v.t.up * S  
    P.t.up <- MortSmooth_BWB(RTBa, RTCt, nba, nbt, v.t.up)
    P.t.up <- kappas[2] * P.t.up
    p.t.up <- MortSmooth_BcoefB(t(Ba), t(Ct), v.t.up*G.t.up)
    p.t.up <- kappas[2] * p.t.up
    ## low
    v.t.low <- BaCt.ALPHAS < G.t.low
    v.t.low <- v.t.low * S  
    P.t.low <- MortSmooth_BWB(RTBa, RTCt, nba, nbt, v.t.low)
    P.t.low <- kappas[2] * P.t.low
    p.t.low <- MortSmooth_BcoefB(t(Ba), t(Ct), v.t.low*G.t.low)
    p.t.low <- kappas[2] * p.t.low
    
    P.t <- P.t.up + P.t.low
    p.t <- p.t.up + p.t.low
    
    ETA <- MortSmooth_BcoefB(Ba, Bt, ALPHAS)
    MU <- E * exp(ETA)
    
    dev.old <- dev
    dev <- 2 * sum(WEI * (Y * log(Y0/MU)))
    d.dev <- abs(dev - dev.old)/(0.1 + abs(dev))
    
    W <- MU
    ZZ <- ETA + (1/MU) * (Y - MU)
    ZZ[which(WEI == 0)] <- 0
    WW <- WEI * W
    tBWB <- MortSmooth_BWB(RTBa, RTBt, nba, nbt, WW)
    tBWBpP <- tBWB + P + P.a + P.t
    tBWZ <- MortSmooth_BcoefB(t(Ba), t(Bt), (WW * ZZ))
    tBWZpP <- tBWZ + p.a + p.t
    alphas <- solve(tBWBpP, c(tBWZpP))
    ALPHAS <- matrix(alphas, nrow = nba)
    if(verbose) cat(it, d.dev, "\n")
    if(d.dev < 10^-4 & it>=4) break   
  }
  ## fitted values
  ALPHAS.hat <- ALPHAS
  ETA.hat <- ETA
  MU.hat <- MU
  ## fitted derivatives over ages
  ETA.hat1a <- MortSmooth_BcoefB(Ca, Bt, ALPHAS.hat)
  ## fitted derivatives over years
  ETA.hat1t <- MortSmooth_BcoefB(Ba, Ct, ALPHAS.hat)
  ## diagnostics
  H <- solve(tBWBpP, tBWB)
  h <- diag(H)
  ed <- sum(h)
  psi <- dev / (sum(WEI) - ed)
  bic <- dev + log(sum(WEI))*ed
  
  
  ## variance-covariance matrix for the coefficients alpha
  tBWB <- MortSmooth_BWB(RTBa, RTBt, nba, nbt, WEI*MU.hat)
  tBWBpP1 <- solve(tBWB + P + P.a/kappas[1] + P.t/kappas[2])# 
  Valphas <- tBWBpP1 %*% tBWB %*% tBWBpP1
  Valphaspsi <- psi*Valphas
  SE.ETA.hat <- Mort2Dsmooth_se(RTBa, RTBt, nba, nbt, BWB.P1=tBWBpP1)
  
  ## return objects
  out <- list(
    ## original data
    Y=Y, E=E, lambdas=lambdas,
    kappas=kappas, deltas=deltas, S=S,
    ## diagnostic
    dev=dev, ed=ed, bic=bic, psi=psi,
    ## fitted values
    ALPHAS=ALPHAS.hat, ETA=ETA.hat, MU=MU.hat,
    ETA1a=ETA.hat1a, ETA1t=ETA.hat1t, 
    ## basis and associated objects
    Ba=Ba, Bt=Bt,
    Ca=Ca, Ct=Ct,
    ## convergence
    d.dev=d.dev,
    it=it,
    ## uncertainty
    Valphas=Valphas, Valphaspsi=Valphaspsi, SE.ETA.hat=SE.ETA.hat
  )
  return(out)
}

## WRAPPER FUNCTION
CCPsplines_fun <- function(ages, cohorts, Z, E, sex,
                           simulate=FALSE, n.sim=100,
                           levels=c(95,50),
                           print.last=TRUE, print.all=FALSE,
                           print.sim=FALSE,
                           yrs.exclude=NULL, decrease=TRUE){

  ## printing the model being fitted
  cat("fitting the CCP-splines model", '\n')
  
  WEI <- matrix(1,nrow(Z),ncol(Z))
  WEI[E==0] <- 0
  WEI[is.na(Z)] <- 0
  ## eventually add weights=0 for user-defined years (e.g. wars)
  if(!is.null(yrs.exclude)){
    # WEI1 <- WEI
    ## matrix with years
    yrs.mat <- row(Z)+col(Z)+c[1]-2
    ## where to exclude
    whi.exclude <- which(yrs.mat%in%yrs.exclude)
    ## include in the weights
    WEI[whi.exclude] <- 0
  }
  ## to get good starting values for the lambdas-search
  if (print.all) cat("Getting starting values ...", "\n")
  FIT0 <- PSinfant(Y=Z, E=E, WEI=WEI, lambdas=c(1,100), verbose=FALSE)
  ## function to extract the BIC for a given lambda
  BIC1 <- function(par){
    FIT1 <- PSinfant(Y=Z, E=E, lambdas=par, WEI=WEI, ALPHAS.st=FIT0$ALPHAS)
    FIT1$bic
  }
  ## optimizing lambdas using greedy grid search
  if (print.all) cat("Optimizing smoothing parameters ...", "\n")
  OPT1 <- cleversearch(BIC1, lower=c(-4, 1), upper=c(0, 5),
                       ngrid=5, logscale=TRUE, verbose=FALSE)
  ## optimal smoothing parameters lambdas
  lambdas.hat <- OPT1$par
  if (print.all) cat("Optimal smoothing parameters:", lambdas.hat, "\n")
  ## estimating mortality with optimal lambdas
  if (print.all) cat("Fitting observed data and computing derivatives ...", "\n")
  FIT1 <- PSinfant(Y=Z, E=E, WEI = WEI, lambdas=lambdas.hat, 
                   ALPHAS.st=FIT0$ALPHAS, verbose=FALSE)
  ## extract deltas
  ETA1a <- FIT1$ETA1a
  ETA1t <- FIT1$ETA1t
  ## discard derivatives in forecast areas
  ETA1a[WEI==0] <- NA
  ETA1t[WEI==0] <- NA
  ## compute levels and deltas over ages
  p.a.up <- (100 - (100-levels[1])/2)/100
  p.a.low <- ((100-levels[1])/2)/100
  delta.a.up <- apply(ETA1a, 1, quantile,
                      probs=p.a.up, na.rm=TRUE)
  delta.a.low <- apply(ETA1a, 1, quantile,
                       probs=p.a.low, na.rm=TRUE)
  ## compute levels and deltas over years
  p.t.up <- (100 - (100-levels[2])/2)/100
  p.t.low <- ((100-levels[2])/2)/100
  delta.t.up <- apply(ETA1t, 1, quantile,
                      probs=p.t.up, na.rm=TRUE)
  delta.t.low <- apply(ETA1t, 1, quantile,
                       probs=p.t.low, na.rm=TRUE)
  if(decrease){
    delta.t.up[delta.t.up>0] <- 0
    delta.t.low[delta.t.low>0] <- 0
  }
  deltas <- list(delta.a.up=delta.a.up,delta.a.low=delta.a.low,
                 delta.t.up=delta.t.up,delta.t.low=delta.t.low)

  
  S <- WEI
  S <- 1-S
  if (print.all) cat("Fitting and forecasting data ...", "\n")
  FIT <- CPSfunction(Y=Z, E=E, WEI=WEI, lambdas=lambdas.hat,
                     deltas=deltas, S=S, verbose=FALSE, ALPHAS.st=FIT0$ALPHAS)
  
  if (print.last) cat("converged at iter.", FIT$it, ", conv. criteria", FIT$d.dev, '\n')
  ## compute e0 and e-dagger for the point estimates
  if (print.all) cat("Computing estimated ex and e-dagger ...", "\n")
  source("funs/LifetableMX.R")
  LE.hat <- ED.hat <- matrix(NA,nrow(Z),ncol(Z)) 
  for (i in 1:ncol(LE.hat)){
    lt <- lifetable.mx(x=x, mx=exp(FIT$ETA)[,i], sex=sex)
    LE.hat[,i] <- lt$ex
    ED.hat[,i] <- eDagger(lt)  
  }
  FIT$LE.hat <- LE.hat
  FIT$ED.hat <- ED.hat
  
  ## adding simulation if required
  if(simulate){
    if (print.all) cat("Simulating coefficients ...", "\n")
    require(MASS)
    ALPHAS.sim <- mvrnorm(n=n.sim, 
                          mu=c(FIT$ALPHAS),
                          Sigma=FIT$Valphas)
    Ba <- FIT$Ba
    Bt <- FIT$Bt
    if (print.all) cat("Computing simulated log-mortality, ex and e-dagger ...", "\n")
    ## empty matrices/arrays to store bootstrap results
    ETA.sim <- LE.sim <- ED.sim <- array(NA,dim=c(nrow(Z),ncol(Z),n.sim))
    ## for each simulation
    for(sim in 1:n.sim){
      if (print.sim & sim %in% round(seq(0,n.sim,length.out = 11))) cat("simulation",sim,"/",n.sim,"\n")
      alphas.s <- ALPHAS.sim[sim,]
      ALPHAS.s <- matrix(alphas.s, ncol(Ba), ncol(Bt))
      ETA.s <- MortSmooth_BcoefB(Ba, Bt, ALPHAS.s)
      ## saving outcomes
      ETA.sim[,,sim] <- ETA.s
      for(i in 1:ncol(Z)){
        LTsim.is <- lifetable.mx(x=x, mx=exp(ETA.s[,i]), sex=sex)
        LE.sim[,i,sim] <- LTsim.is$ex
        ED.sim[,i,sim] <- eDagger(LTsim.is)
      }
    }
    sim <- list(LE.sim=LE.sim, ED.sim=ED.sim, ETA.sim=ETA.sim)
    out <- append(FIT,sim)
    return(out)
  }else{
    ## output
    return(FIT)
  }
}
