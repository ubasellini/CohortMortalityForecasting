
## --------------------------------------------------------- ##
##
## Functions for fitting the Cohort LLC model
##
## --------------------------------------------------------- ##


## wrapper function for fitting cohort LLC model 
## (allowing to obtain simulations for CIs)
cohort_LLC_fun <- function(ages,cohorts,Z,E,sex,
                           Alpha=NULL,Beta=NULL,KappaL=NULL,
                           simulate=F,n.sim=100,
                           max.iter=1000,tol=1e-05,
                           print.last=T,print.all=F,print.sim=F){
  
  ## printing the model being fitted
  cat("fitting the cohort LLC model", '\n')
  
  ## start by fitting cohort LLC model
  fit <- fit_cohort_LLC_fun(ages = ages,cohorts = cohorts,Z=Z,E=E,sex=sex,
                            Alpha = Alpha,Beta=Beta,KappaL=KappaL,
                            max.iter = max.iter,tol=tol,
                            print.last = print.last,print.all = print.all)
  
  ## adding simulation if required
  if (simulate){
    ## dimensions of interest
    m <- length(ages)
    nc <- length(cohorts)
    
    ## empty matrices/arrays to store simulation results
    Alpha.sim <- Beta.sim <- matrix(NA,m,n.sim)
    Kappa.sim <- matrix(NA,nc,n.sim)
    ETA.sim <- LE.sim <- ED.sim <- Z.sim <- array(NA,dim=c(m,nc,n.sim))
    
    sim <- 1
    for (sim in 1:n.sim){
      ## printing simulation number if required
      if (print.sim & sim %in% round(seq(0,n.sim,length.out = 11))) cat("simulation",sim,"/",n.sim,"\n")
      
      ## simulate Poisson deaths from observed deaths
      Zsim <- suppressWarnings(matrix(rpois(n=m*nc,lambda=c(Z)),m,nc))
      # image.plot(c,x,t(cZ))
      # image.plot(c,x,t(cZsim))
      ## refitting LC
      fit_cohort_LLC_sim <- fit_cohort_LLC_fun(ages=ages,cohorts=cohorts,Z=Zsim,E=E,
                                               sex=sex,print.last = F)
      ## saving outcomes
      Alpha.sim[,sim] <- fit_cohort_LLC_sim$Alpha
      Beta.sim[,sim] <- fit_cohort_LLC_sim$Beta
      Kappa.sim[,sim] <- fit_cohort_LLC_sim$Kappa
      LE.sim[,,sim] <- fit_cohort_LLC_sim$LE.hat
      ED.sim[,,sim] <- fit_cohort_LLC_sim$ED.hat
      ETA.sim[,,sim] <- fit_cohort_LLC_sim$ETA
      Z.sim[,,sim] <- fit_cohort_LLC_sim$Z.hat
    }
    
    ## output
    sim <- list(Alpha.sim=Alpha.sim,Beta.sim=Beta.sim,Kappa.sim=Kappa.sim,
                LE.sim=LE.sim,ED.sim=ED.sim,ETA.sim=ETA.sim,Z.sim=Z.sim)
    out <- append(fit,sim)
    return(out)
    
  }else{
    ## output
    return(fit)
  }
  
}




## Starting values of LC model
LC_starting_pars <- function(Z,E){
  ## dimensions
  m <- nrow(Z)
  n <- ncol(Z)
  ## starting values
  One <- matrix(1, nrow=n, ncol = 1)    
  Fit <- log((Z + 1)/(E + 2))
  ## for Alpha, take mean of log death rates
  Alpha <- apply(Fit/n,1,sum,na.rm=T)
  ## for Beta, take alpha and normalize to sum to one
  Beta <- matrix(1 * Alpha, ncol = 1)
  sum.Beta <- sum(Beta) 
  Beta <- Beta / sum.Beta
  ## for Kappa, standardize a series from n to 1 
  Kappa <- matrix(seq(n, 1, by = -1), nrow = n, ncol = 1)
  Kappa <- Kappa - mean(Kappa)
  Kappa <- Kappa / sqrt(sum(Kappa * Kappa))
  ## return list
  out <- list(Alpha=Alpha,Beta=Beta,Kappa=Kappa,One=One)
}


## function for fitting the cohort LC model 
fit_cohort_LLC_fun <- function(ages,cohorts,Z,E,sex,
                               Alpha=NULL,Beta=NULL,KappaL=NULL,
                               max.iter=1000,tol=1e-05,
                               print.last=T,print.all=F){
  ## dimensions
  m <- nrow(Z)
  n <- ncol(Z)
  mn <- m*n
  ## time points
  t <- 1:n
  ## ages, cohorts and periods matrices
  AGES <- matrix(ages,nrow=m,ncol=n)
  COHORTS <- matrix(rep(cohorts,each=m),nrow=m)
  PERIODS <- COHORTS + AGES
  ## derive periods
  periods <- min(PERIODS):max(PERIODS)
  last.per.obs <- cohorts[n]+sum(!is.na(Z[,n]))-1
  periods.obs <- periods[1]:last.per.obs
  ng <- length(periods)
  ng.obs <- length(periods.obs)
  ## Weights: no exposures, missing data, dummy periods
  W <- matrix(1,m,n)
  W[E==0] <- 0
  W[is.na(E)] <- 0
  ## get starting values
  fict.Z <- 10^3
  One <- matrix(1, nrow=n, ncol = 1)
  if (is.null(Alpha)){
    StartingPars <- LC_starting_pars(Z,E)
    Alpha <- StartingPars$Alpha
    Beta <- StartingPars$Beta
  }
  KappaL <- c(1,-0.02)
  ## Model matrix and constraint for Kappa (Linear model)
  X <- cbind(rep(1,n),t)
  KappaL[1] <- -KappaL[2]*sum(t)/n
  ## plotting
  # plot(x,Alpha);plot(x,Beta);plot(c,Kappa)
  
  ETA <- Alpha %*% t(One) + Beta %*% t(X%*%KappaL)
  Z.hat <- E * exp(ETA)
  Z.hat[is.na(Z.hat)] <- fict.Z

  ## fitting model
  fit <- cohort_LLC_iterations(x=ages,Alpha=Alpha,Beta=Beta,KappaL=KappaL,ETA=ETA,W=W,
                               Z.hat=Z.hat,Z=Z,E=E,sex=sex,
                               max.iter=max.iter,tol=tol,print.last=print.last,print.all=print.all)
  ## output
  return(fit)
}

cohort_LLC_iterations <- function(x,Alpha,Beta,KappaL,ETA,W,Z.hat,Z,E,sex,
                                 max.iter,tol,
                                 print.last,print.all){
  
  ## dimensions
  m <- nrow(Z)
  n <- ncol(Z)
  ## time points
  t <- 1:n
  ## Model matrix 
  X <- cbind(rep(1,n),t)
  ## fictitious deaths and vector One
  One <- matrix(1, nrow=n, ncol = 1)
  fict.Z <- 10^3
  ## convergence
  conv <- T
  ## iterations
  for (iter in 1:max.iter){
    Alpha.old <- Alpha
    Beta.old <- Beta
    KappaL.old <- KappaL
    ETA.old <- ETA
    
    ## update alphas 
    c11 <- apply(W*Z.hat, 1, sum) ## sum of fitted deaths for each age
    ## model matrix for alphas
    Ca <- diag(as.vector(c11))
    Z.diff <- Z - Z.hat
    Z.diff[is.na(Z.diff)] <- fict.Z
    ra <- apply(W*Z.diff, 1, sum) ## sum of the difference between actual and fitted deaths for each age
    tryAlpha <- try(solve(Ca, ra + Ca %*% Alpha.old),silent=T)
    if(class(tryAlpha)[1]=="try-error" | any(is.na(tryAlpha)) | all(tryAlpha==0)){
      conv <- F
      break
    }else{
      Alpha <- c(tryAlpha)
    }
    
    ## fitted linear predictor 
    ETA <- Alpha %*% t(One) + Beta %*% t(X%*%KappaL)
    ## fitted expected values 
    Z.hat <- E * exp(ETA)
    Z.hat[is.na(Z.hat)] <- fict.Z
    
    ## update betas
    c22 <- (W*Z.hat) %*% (X%*%KappaL)^2 ## actual deaths by kappa-squared
    ## model matrix for betas
    Cb <- diag(as.vector(c22))
    Z.diff <- Z - Z.hat
    Z.diff[is.na(Z.diff)] <- fict.Z
    rb <- (W*Z.diff) %*% (X%*%KappaL) ## sum of the (difference between actual and fitted deaths multiplied by kappa) for each age
    tryBeta <- try(solve(Cb, rb + Cb %*% Beta.old),silent=T)
    if(class(tryBeta)[1]=="try-error" | any(is.na(tryBeta)) | all(tryBeta==0)){
      conv <- F
      break
    }else{
      Beta <- c(tryBeta)
    }
    
    ## adding the constraints
    sum.Beta <- sum(Beta) 
    Beta <- Beta / sum.Beta
    
    ## fitted linear predictor 
    ETA <- Alpha %*% t(One) + Beta %*% t(X%*%KappaL)
    ## fitted expected values 
    Z.hat <- E * exp(ETA)
    Z.hat[is.na(Z.hat)] <- fict.Z
    
    ## plotting
    # plot(Beta.old,ylim=range(Beta.old,Beta));points(Beta,col=2)
    
    ## update Kappa with Newton-Raphson
    Ck <- t(X) %*% diag(as.vector(t(W*Z.hat)%*%(Beta^2))) %*% X 
    Z.diff <- Z - Z.hat
    Z.diff[is.na(Z.diff)] <- fict.Z
    rk  <- t(X) %*% t(W*Z.diff) %*% Beta  ## ## sum of the (difference between actual and fitted deaths multiplied by beta) for each year
    tryKappaL   <- try(solve(Ck, rk + Ck %*% KappaL.old),silent=T)
    if(class(tryKappaL)[1]=="try-error" | any(is.na(tryKappaL)) | all(tryKappaL==0)){
      conv <- FALSE
      break
    }else{
      KappaL <- c(tryKappaL)
    }
    
    ## constraint for kappas
    KappaL[1] <- -KappaL[2]*sum(t)/n
    
    ## fitted linear predictor 
    ETA <- Alpha %*% t(One) + Beta %*% t(X%*%KappaL)
    ## fitted expected values 
    Z.hat <- E * exp(ETA)
    Z.hat[is.na(Z.hat)] <- fict.Z
    
    ## plotting 
    # plot(c(t(X%*%KappaL.old)),ylim=range(t(X%*%KappaL.old),t(X%*%KappaL)));points(c(t(X%*%KappaL)),col=2)
    
    ## break if reaching max iter
    if (iter==max.iter){
      conv <- F
      break
    }
    
    ## torelance criterion
    crit <- max(abs(ETA - ETA.old))
    if (crit < tol & conv){
      if (print.last) cat("converged at iter.",iter,", conv. criteria", crit, '\n')
      break
    } 
    if (print.all) cat(iter, crit, '\n')
  }
  
  if (conv){
    
    ## update linear predictor and fitted deaths
    Kappa <- t(X%*%KappaL)
    ## fitted linear predictor 
    ETA <- Alpha %*% t(One) + Beta %*% t(c(Kappa))
    ## fitted expected values 
    Z.hat <- E * exp(ETA)
    
    ## compute life exp and e-dagger
    ## taking observed rates where possible and forecast ones
    source("funs/LifetableMX.R")
    ## forecasting region
    FORE <- matrix(1,m,n)
    FORE[W==1] <- 0
    MXLC <- Z/E
    MXLC[FORE==1] <- exp(ETA[FORE==1])
    LE.hat <- ED.hat <- matrix(NA,m,n) 
    for (i in 1:n){
      lt <- lifetable.mx(x=x,mx=MXLC[,i],sex=sex)
      LE.hat[,i] <- lt$ex
      ED.hat[,i] <- eDagger(lt)  
    }
    ## compute deviance, ED, BIC 
    DEV <- 2 * sum(W * Z * log(ifelse(Z==0, 1e-08, Z) / ifelse(Z.hat==0, 1e-08, Z.hat)),na.rm = T)
    EDpars <- 2*m + 2
    ED <- EDpars - 2   ## remove two constrains 
    BIC <- DEV + log(sum(W)) * ED
    ## compute deviance residuals
    t1 <- sign(Z-Z.hat)
    t2 <- 2*(W * Z * log(ifelse(Z==0, 1e-08, Z) / ifelse(Z.hat==0, 1e-08, Z.hat))-W*(Z-Z.hat))
    DEVres <- t1*sqrt(t2)
    ## output
    out <- list(Alpha=Alpha,Beta=Beta,Kappa=Kappa,ETA=ETA,Z.hat=Z.hat,LE.hat=LE.hat,
                ED.hat=ED.hat,
                DEV=DEV,ED=ED,BIC=BIC,DEVres=DEVres,iter=iter,crit=crit,conv=conv)
  }else{
    if (print.last) cat("no convergence", '\n')
    out <- list(conv=conv)
  }
  return(out)
}



