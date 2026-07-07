


## --------------------------------------------------------- ##
##
## Functions for fitting the diagonal LC model
##
## --------------------------------------------------------- ##



## wrapper function to obtain cohort forecasts 
## from a period LC forecasts
# ages=x
# years=t1
# cohorts=c1
# Z=Z1
# E=E1
# cZ=cZ1
# cE=cE1
# sex=sex
# n.sim=n.sim
# Alpha=NULL
# Beta=NULL
# Kappa=NULL
# lambdaA=0
# lambdaB=0
# max.iter=1000
# tol=1e-05
# print.last=T
# print.all=F


diagonal_LC_fun <- function(ages,years,cohorts,
                            Z,E,cZ,cE,sex,
                            Alpha=NULL,Beta=NULL,Kappa=NULL,
                            lambdaA=0,lambdaB=0,
                            n.sim=100,
                            max.iter=1000,tol=1e-05,
                            print.last=T,print.all=F){
  
  
  ## printing the model being fitted
  cat("fitting the diagonal LC model", '\n')
  
  ## dimensions
  m <- length(ages)
  n <- length(years)
  nc <- length(cohorts)
  ## find number of forecast years needed to complete cohorts
  t.fore <- (years[n]+1):(cohorts[nc]+m+10)
  n.fore <- length(t.fore)
  
  ## fit period LC
  fit <- fit_period_LC_fun(ages = ages,years = years,
                           Z=Z,E=E,sex=sex,n.fore=n.fore,
                           n.sim=n.sim,
                           Alpha = Alpha,Beta=Beta,Kappa=Kappa,
                           max.iter = max.iter,tol=tol,
                           print.last = print.last,print.all = print.all)
  
  ## extract cohort forecasts
  # library(fields)
  ## observed rates
  MX <- Z/E
  cMX <- cZ/cE
  
  ## fitted rates
  MX.fit <- exp(fit$ETA)
  # image.plot(t,x,t(log(MX)))
  # image.plot(t,x,t(log(MX.fit)))
  cMX.fit <- matrix(NA,nrow(cMX),ncol(cMX))
  ## for the last observed year, cohort data can't be computed as it needs 
  ## the first forecast. Hence include the median forecast for the first forecast year
  MX.fore1 <- exp(apply(fit$ETA.sim[,1,],1,median))
  MX.fit.ext <- cbind(MX.fit,MX.fore1)
  # image(c(t,2020),x,t(log(MX.fit.ext)))
  AGES <- matrix(ages,nrow=m,ncol=n+1)
  PERIODS <- matrix(rep(c(years,years[n]+1),each=m),nrow=m)
  COHORTS  <-  PERIODS - AGES
  c.temp <- COHORTS[1,1]:max(cohorts)
  nc.temp <- length(c.temp)
  j <- i <- 1
  for (i in 1:nc.temp){
    my.coh <- which(c==c.temp[i])
    for (j in 1:m){
      whi <- which(COHORTS[j,]==c.temp[i])
      whi1 <- which(COHORTS[j,]==c.temp[i]+1) 
      if (length(whi)!=0 & length(whi1)!=0){
        mx.temp <- MX.fit.ext[j,whi]
        mx.temp1 <- MX.fit.ext[j,whi1]
        ## using Carl's adjustment (AP2 estimator)
        cMX.fit[j,my.coh] <- (mx.temp + mx.temp1) / 2
      }
    }
  }
  # image.plot(c,x,t(log(cMX)))
  # image.plot(c,x,t(log(cMX.fit)))
  
  ## all time periods, and APC matrices
  t <- years
  t.all <- c(t,t.fore)
  nt.all <- length(t.all)
  AGES <- matrix(ages,nrow=m,ncol=nt.all)
  PERIODS <- matrix(rep(t.all,each=m),nrow=m)
  COHORTS  <-  PERIODS - AGES
  
  ## objects to store results
  ETA <- LE <- ED <- array(NA,dim=c(m,nc,n.sim))
  sim <- 1
  for (sim in 1:n.sim){
    MX.temp <- cbind(MX.fit,exp(fit$ETA.sim[,,sim]))
    cMX.temp <- matrix(NA,m,nc)
    i <- 1
    for (i in 1:nc){
      whi <- which(COHORTS==c[i])
      whi1 <- which(COHORTS==c[i]+1) 
      mx.temp <- MX.temp[whi]
      mx.temp1 <- MX.temp[whi1]
      if (length(mx.temp)<m){
        mx.temp <- c(rep(NA,m-length(mx.temp)),mx.temp)
        mx.temp1 <- c(rep(NA,m-length(mx.temp1)),mx.temp1)
      } 
      ## using Carl's adjustment (AP2 estimator)
      cMX.temp[,i] <- (mx.temp + mx.temp1) / 2
    }
    
    ## create dataframe to store forecast and select only forecast
    ## (i.e. keep observed data otherwise)
    df.temp <- tibble(ages=rep(ages,nc),cohorts=rep(cohorts,each=m),
                      observed=c(cMX.fit),forecast=c(cMX.temp)) %>% 
      mutate(rates=case_when(
        is.na(observed) ~ forecast,
        TRUE            ~ observed
      ))
    
    ETA[,,sim] <- log(df.temp$rates)

    ## life-table
    for (i in 1:nc){
      if (!is.na(ETA[1,i,sim])){
        clt <- lifetable.mx(x=x,mx=exp(ETA[,i,sim]),sex=sex)
        LE[,i,sim] <- clt$ex
        ED[,i,sim] <- eDagger(clt)
      }
    }
    
  }
  
  ## output
  out <- list(LE.sim=LE,ED.sim=ED,ETA.sim=ETA)
  return(out)
  
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
fit_period_LC_fun <- function(ages,years,
                              Z,E,sex,n.sim=100,n.fore,
                              lambdaA=0,lambdaB=0,
                              Alpha=NULL,Beta=NULL,Kappa=NULL,
                              max.iter=1000,tol=1e-05,
                              print.last=T,print.all=F){
  ## dimensions
  m <- nrow(Z)
  n <- ncol(Z)
  mn <- m*n
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
    Kappa <- StartingPars$Kappa
  }
  ETA <- Alpha %*% t(One) + Beta %*% t(Kappa)
  Z.hat <- E * exp(ETA)
  Z.hat[is.na(Z.hat)] <- fict.Z

  ## fitting model
  fit <- period_LC_iterations(ages=ages,years=years,n.sim=n.sim,n.fore=n.fore,
                              Alpha=Alpha,Beta=Beta,Kappa=Kappa,ETA=ETA,W=W,
                              Z.hat=Z.hat,Z=Z,E=E,sex=sex,
                              lambdaA=lambdaA,lambdaB=lambdaB,lambdaK=0,
                              max.iter=max.iter,tol=tol,print.last=print.last,print.all=print.all)
  
  return(fit)
}

period_LC_iterations <- function(ages,years,n.sim,n.fore,
                                 Alpha,Beta,Kappa,ETA,W,Z.hat,Z,E,sex,
                                 lambdaA,lambdaB,lambdaK,
                                 max.iter,tol,
                                 print.last,print.all){
  
  ## dimensions
  m <- nrow(Z)
  n <- ncol(Z)
  ## penalties and difference-matrices (for smooth LC)
  Dm <- diff(diag(m), diff=2)
  Dn <- diff(diag(n), diff=2)
  Pa <- lambdaA * (t(Dm) %*% Dm) 
  Pb <- lambdaB * (t(Dm) %*% Dm)
  Pk <- lambdaK * (t(Dn) %*% Dn)
  ## fictitious deaths and vector One
  One <- matrix(1, nrow=n, ncol = 1)
  fict.Z <- 10^3
  ## convergence
  conv <- T
  ## iterations
  for (iter in 1:max.iter){
    Alpha.old <- Alpha
    Beta.old <- Beta
    Kappa.old <- Kappa
    ETA.old <- ETA
    
    ## update alphas 
    c11 <- apply(W*Z.hat, 1, sum) ## sum of fitted deaths for each age
    ## model matrix for alphas
    Ca <- diag(as.vector(c11))
    Z.diff <- Z - Z.hat
    Z.diff[is.na(Z.diff)] <- fict.Z
    ra <- apply(W*Z.diff, 1, sum) ## sum of the difference between actual and fitted deaths for each age
    tryAlpha <- try(solve(Ca + Pa, ra + Ca %*% Alpha.old),silent=T)
    if(class(tryAlpha)[1]=="try-error" | any(is.na(tryAlpha)) | all(tryAlpha==0)){
      conv <- F
      break
    }else{
      Alpha <- c(tryAlpha)
    }
    
    ## fitted linear predictor 
    ETA <- Alpha %*% t(One) + Beta %*% t(Kappa)
    ## fitted expected values 
    Z.hat <- E * exp(ETA)
    Z.hat[is.na(Z.hat)] <- fict.Z
    
    ## update betas
    c22 <- (W*Z.hat) %*% Kappa^2 ## actual deaths by kappa-squared
    ## model matrix for betas
    Cb <- diag(as.vector(c22))
    Z.diff <- Z - Z.hat
    Z.diff[is.na(Z.diff)] <- fict.Z
    rb <- (W*Z.diff) %*% Kappa ## sum of the (difference between actual and fitted deaths multiplied by kappa) for each age
    tryBeta <- try(solve(Cb + Pb, rb + Cb %*% Beta.old),silent=T)
    if(class(tryBeta)[1]=="try-error" | any(is.na(tryBeta)) | all(tryBeta==0)){
      conv <- F
      break
    }else{
      Beta <- c(tryBeta)
    }
    
    ## fitted linear predictor 
    ETA <- Alpha %*% t(One) + Beta %*% t(Kappa)
    ## fitted expected values 
    Z.hat <- E * exp(ETA)
    Z.hat[is.na(Z.hat)] <- fict.Z
    
    ## model matrix for kappas
    Ck <- diag(as.vector(t(W*Z.hat)%*%(Beta^2)))
    Z.diff <- Z - Z.hat
    Z.diff[is.na(Z.diff)] <- fict.Z
    rk  <- t(W*Z.diff) %*% Beta ## ## sum of the (difference between actual and fitted deaths multiplied by beta) for each year
    tryKappa   <- try(solve(Ck + Pk, rk + Ck %*% Kappa.old),silent=T)
    if(class(tryKappa)[1]=="try-error" | any(is.na(tryKappa)) | all(tryKappa==0)){
      conv <- F
      break
    }else{
      Kappa <- c(tryKappa)
    }
    
    ## constraint for kappas
    Kappa <- Kappa - mean(Kappa)
    Kappa <- Kappa / sqrt(sum(Kappa*Kappa))
    
    ## fitted linear predictor 
    ETA <- Alpha %*% t(One) + Beta %*% t(Kappa)
    ## fitted expected values 
    Z.hat <- E * exp(ETA)
    Z.hat[is.na(Z.hat)] <- fict.Z
    
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
    ## constraints
    sum.Beta <- sum(Beta)
    Beta <- Beta / sum.Beta
    Kappa <- Kappa * sum.Beta
    ## fitted linear predictor 
    ETA <- Alpha %*% t(One) + Beta %*% t(Kappa)
    ## fitted expected values 
    Z.hat <- E * exp(ETA)
    
    ## compute life exp and e-dagger
    ## taking observed rates where possible and forecast ones
    source("funs/LifetableMX.R")
    ## forecasting region
    MXLC <- exp(ETA)
    LE.hat <- ED.hat <- matrix(NA,m,n) 
    for (i in 1:n){
      lt <- lifetable.mx(x=ages,mx=MXLC[,i],sex=sex)
      LE.hat[,i] <- lt$ex
      ED.hat[,i] <- eDagger(lt)  
    }
    ## compute deviance, ED, BIC 
    DEV <- 2 * sum(W * Z * log(ifelse(Z==0, 1e-08, Z) / ifelse(Z.hat==0, 1e-08, Z.hat)),na.rm = T)
    EDpars <- sum(diag(solve(Ca + Pa,Ca))) + sum(diag(solve(Cb + Pb, Cb))) + sum(diag(solve(Ck+ Pk, Ck)))
    ED <- EDpars - 2   ## remove two constrains 
    BIC <- DEV + log(sum(W)) * ED
    ## compute deviance residuals
    t1 <- sign(Z-Z.hat)
    t2 <- 2*(W * Z * log(ifelse(Z==0, 1e-08, Z) / ifelse(Z.hat==0, 1e-08, Z.hat))-W*(Z-Z.hat))
    DEVres <- t1*sqrt(t2)
    
    ##--- forecasting with LC simulations
    require(forecast)
    RWD <- Arima(ts(Kappa,start = years[1]), order=c(0,1,0), include.drift=TRUE) 
    
    ## define matrices to store simulation results
    One.fore <- matrix(1, nrow=n.fore, ncol = 1)
    Kappa.sim <- matrix(NA,nrow=n.fore,ncol=n.sim)
    ETA.sim <- array(NA,dim=c(m,n.fore,n.sim))
    
    ## simulate future K
    i <- 1
    for(i in 1:n.sim){
      ## generate simulation with bootsrapping
      kappa.sim <- simulate(RWD, nsim=n.fore,
                            future=TRUE, bootstrap=TRUE)
      # plot(y1,kappa,ylim=range(kappa,kappa.sim),xlim=range(y))
      # lines(yF,kappa.sim)
      
      ## derive the bootsrap values
      Kappa.sim[,i] <- kappa.sim
      ETA.sim[,,i] <- Alpha %*% t(One.fore) + Beta %*% t(kappa.sim)
    }
    
    ## output
    out <- list(Alpha=Alpha,Beta=Beta,Kappa=Kappa,ETA=ETA,Z.hat=Z.hat,LE.hat=LE.hat,
                ED.hat=ED.hat,lambda=lambdaA,Kappa.sim=Kappa.sim,ETA.sim=ETA.sim,
                DEV=DEV,ED=ED,BIC=BIC,DEVres=DEVres,iter=iter,crit=crit,conv=conv)
  }else{
    if (print.last) cat("no convergence", '\n')
    out <- list(conv=conv)
  }
  return(out)
}



