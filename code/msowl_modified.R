# Multistate outcome weighted learning
# Bakoyannis Biometrika 2023 

# libraries 
library(survival)
library(kernlab)
# library(pbapply)  # parallel
# library(parallel)
library(future)
library(future.apply)
library(matrixcalc)
library(dplyr)
#setwd("H:/msowl_hiv/programs")
#setwd("/Users/giorgosbakoyannis/Documents/MSOWL_HIV")
#dat <- read.csv("example_data.csv")

# plan(multisession)  # Use multisession for parallel processing on the local machine
# lambda.vec = c(0.00001, 0.0001, 0.001, 0.01, 0.1, 0.5, 1, 5, 10, 50, 100, 1000, 10000)
lambda.vec = c(0.001, 0.01, 0.1, 0.5, 1, 5, 10, 50, 100)
# Weighted SVM solver
# Code adapted from DTRlearn2 package (source code:
# https://github.com/ychen178/DTRlearn2/blob/master/R/wsvm_solve.R)
wsvm_solve <-function(X, A, wR, kernel='linear', sigma=0.05, lambda=1, e=1e-7) {
  if (kernel=='linear') {
    K = X %*% t(X)
    if (is.vector(X)) K = t(X) %*% X
  }
  else if (kernel=='rbf'){
    rbf = rbfdot(sigma = sigma)
    K = kernelMatrix(rbf, X)
  }
  y = A * sign(wR)
  H = y %*% t(y) * K
  n = length(A)
  C = 1/(2*n*lambda)
  solution <-
  tryCatch(
    ipop(
       c = rep(-1, n), 
       H = H, 
       A = t(y), 
       b = 0, 
       l = numeric(n), 
       u = C*abs(wR), 
       r = 0
      ), 
  error=function(er) er
  )
  
  if ("error" %in% class(solution)) {
    model <- NA
    class(model) <- "wsvm_solve:ipop function error"
  } else {
    alpha = primal(solution)
    alpha1 = alpha * y

    if (kernel=='linear'){
      w = t(X) %*% alpha1
      fitted = X %*% w
    } else if (kernel=='rbf'){
      fitted = K %*% alpha1
    }
    rm = y - fitted
    Imid = (alpha < C-e) & (alpha > e)
    rmid = rm[Imid==1]
    if (sum(Imid)>0){
      bias = mean(rmid)
    } else {
      Iup = ((alpha<e)&(A==-sign(wR)))|((alpha>C-e)&(A==sign(wR)))
      Ilow = ((alpha<e)&(A==sign(wR)))|((alpha>C-e)&(A==-sign(wR)))
      rup = rm[Iup]
      rlow = rm[Ilow]
      bias = (min(rup)+max(rlow))/2
    }

    if (kernel=='linear') {
      model = list(beta0=bias, beta=w, alpha1=alpha1)
      class(model)<-'linear'
    } else if (kernel=='rbf') {
      model = list(beta0=bias, sigma=sigma, Z=X, alpha1=alpha1)
      class(model) = 'rbf'
    }
  }
  return (model)
}


# zeta_t() --------
# Utility function for the optimal treatment rule estimation function
# Calc the bracketed expr from the Value / Empirical Risk function 
# \int_0^{\tau} { \frac{Y_w(t)I(C >= min(T,T))}{exp(-\Lambda_0(min(\tilde{T},t))} dm(t) } 
# which calcs the Mean Patient-Weighted Multistate Survival Time, 
# Inflated by Inverse Probability of Censoring Survival Times 

zeta_t <- function(t, data, w=c(0, 1, 0), hazard, S, T){
  
  # patient death indicator per state before time t
  events <- list()
  for(j in 1:length(S)){
    if(!(S[j] %in% T)){
      # not terminal state 
      events[[j]] <- (data$s1==S[j] & data$t1<=t & data$t2>t)
    } else {
      # terminal state 
      events[[j]] <- (data$s2==S[j] & data$t2<=t)
    }
  }
  # weight the events by patient-preference of state
  event <- do.call(cbind, events)%*%w

  Tt <- (1*!(data$s2 %in% T))*t + (1*(data$s2 %in% T))*mapply(min, data$t2, t) #T^t

  # inverse probability of censoring weighting  
  # get the element of the instantaneous hazard rate at each time Tt
  element <- function(t){
    max(1,max(1:length(hazard$time)*(hazard$time<=t)))
  }
  elements <- sapply(Tt, element)
  G_t <- exp(-hazard$hazard[elements])
  res <- event/G_t
  
  return(res)
}

ms <- function(t1, t2, s1, s2, a){
  #if (missing(id))
  #  stop("Missing id argument")
  if (missing(t1))
    stop("Missing values in `t1'")
  if (missing(t2))
    stop("Missing values in `t2'")
  if (missing(s1))
    stop("Missing values in `s1'")
  if (missing(s2))
    stop("Missing values in `s2'")
  if (max(t1 >= t2)==1)
    stop("t1 >= t2")
  if (max(!(unique(a) %in% c(-1, 1)))==1)
    stop("Treatment var `a' is not equal to -1 or 1")
  msdat <- cbind(t1 = t1, t2 = t2, s1 = s1, s2 = s2, a = a)
  class(msdat) <- "msdat"
  msdat
}


calculate_reward = function(
    formula, 
    id, 
    w = c(0, 1, 0), 
    tau = 3, 
    data,
    owl_method = "MSOWL"
) {
  msout <- model.frame(formula, data)[[1]]
  covs <- model.matrix(formula, data)[, -1]
  data <- cbind.data.frame(id = data$id, as.data.frame(cbind(msout, covs)))
  feat <- colnames(data)[(ncol(msout) + 2):ncol(data)]
  rm(msout, covs)
  data <- data[order(data[,id], data[,"t1"]),]
  ## state space
  S <- with(data, sort(unique(c(s1, s2))))
  ## absorbing state subspace
  T <- S[!(S %in% sort(unique(data$s1)))]
  # validate state patient-preference weights 
  if(length(S)!=length(w)){
    stop("length(w) not equal to the total number of states")
  }
  if(min(w)<0 | max(w)>1 | sum(w)>=length(w)){
    stop("min(w)<0 or max(w)>1 or sum(w)>=length(w)")
  }
  
  # Calculate reward
  # Censoring Distribution. 
  c.delta <- 1*(data$s1==data$s2 & data$t2<tau)
  form <- as.formula("Surv(time=t1, time2=t2, event=c.delta, type='counting') ~ 1")
  # nelson-aalen estimator for baseline hazard (via coxph model) 
  fit <- coxph(form, data=data)
  hazard <- basehaz(fit, centered=FALSE) # baseline hazard function 
  
  # MSOWL INTEGRATES THE SURVIVAL OF EACH PATIENT OVER TIME UNTIL CENSOR TIME 
  # REQUIRES MORE COMPUTATIONALLY COSTLY ANALYSIS OF HAZARD 
  # USES FUNCTION zeta_t() 
  if(owl_method == "MSOWL") {
    # combine 
    # (1) imputed survival times from both censored and uncensored data
    # (2) true un-censored survival times 
    tt <- sort(
      unique(
        c(
          survfit(fit, se.fit = FALSE)$time, # censoring times
          data[data$s1!=data$s2, "t2"]       # state transition times 
        )
      )
    )
    tt <- c(0,tt[tt<=tau]) # include zero and only within tau 
    
    # zeta_t: compute hazard function (all pairwise combinations )
    event <- matrix(
      # pbsapply(
      sapply(
          tt, 
          zeta_t, 
          data=data, 
          w=w, 
          hazard=hazard, 
          S=S, 
          T=T
        ),
        byrow=TRUE, 
        nrow=length(tt)
      )
    
    # per-individual integration 
    # hazard(t) = events(t) * dm(t)
    dm <- diff(c(tt, tau), lag=1) # dm 
    xi <- colSums(event*dm)  
    # length(xi)
    
    # this just transposes from row to column ... 
    xi <- matrix(tapply(X=xi, INDEX = data$id, FUN=sum))
    
  } else if (owl_method == "ICO") {
    # INVERSE CENSORING OUTCOME (ICO) from ZhaoKosorok2015 Doubly Robust paper 
    # # inverse probability of censoring weighting  
    element <- function(t){
      # # get the element of the instantaneous hazard rate at each time Tt
      max(1,max(1:length(hazard$time)*(hazard$time<=t)))
    }
    elements <- sapply(data$t2, element)
    G_t <- exp(-hazard$hazard[elements])
    xi.ico <- (!c.delta*(data$t2 - data$t1))/G_t
    xi.ico <- matrix(tapply(X=xi.ico, INDEX = data$id, FUN=sum))
    xi = xi.ico
    # lala = setDT(
    #   data.frame(
    #     msowl = xi,
    #     ico = xi.ico,
    #     G_t = G_t,
    #     t2 = data$t2
    #   )
    # )
    # lala[, match := msowl == ico]
    # lala[match == F, ]
    
  } else if (owl_method == "DR") {
    # ZhaoKosorok2015 Doubly Robust paper 
    # DR: W_i = R(Y_i, ind_cens_i, \hat{S_C},\hat{E_T})
    # R(Y,\delta,S_C,E_\tilde{E}) = \frac{\Delta Y}{S_C(Y|A,X)}-\int{E_\tilde{T}}
    # model survival function S_T(T|A,X) 
    # model E_T(T_tilde | T > t, A, X)
    # model S_C(t|A,X)  censoring distribution
    # W_i = R(Y_i, \Delta_i, S_C, E_T)
  }
  
  # propensity score 
  dat <- data[data$t1==0,]
  dat <- dat[order(dat[, id]),]
  pi_n <- mean(dat$a==1)
  
  # outcome weighting
  Wt <- xi/(pi_n*(dat$a==1) + (1 - pi_n)*(dat$a==-1))
  
  complete_data = list(
    id = dat$id,
    X = as.matrix(dat[,feat]),
    A = dat$a, 
    Y = dat$t2 - dat$t1, 
    c.delta = c.delta,
    Wt = Wt
  )
  
  retlist = list(
    dat, 
    complete_data,
    feat
  )
  
  return(retlist)
}


#Estimator of the value function of an ITR
#function(formula, id, w = c(0, 1, 0), tau = 3, data, lambda = 1){
msowl.val <- function(
    formula, 
    id, 
    w = c(0, 1, 0), 
    tau = 3, 
    data,
    rule, 
    fixed.rules = FALSE, 
    SE = FALSE, 
    nboot = 100,
    trim = NULL,
    reward = NULL,
    debug = F
){
  # print(rule)
  # print(class(rule))
  # if( !(class(rule) %in%  c("msowl", "data.frame") ) ) 
    # print(paste(typeof(rule), class(rule)))
    # stop("`rule' should either be output of msowl() run, or a data.frame of {id, treatments}")
  if(!(fixed.rules %in% c(TRUE, FALSE)))
    stop("`fixed.rules' should be TRUE or FALSE")
  if(!(SE %in% c(TRUE, FALSE)))
    stop("`SE' should be TRUE or FALSE")
  if(!is.null(trim)){
    if(trim < 0 | trim > 0.5){
      stop("`trim' should be between >=0 and <= 0.5")
    }
  }
  V_d <- function(data){
  
    if (is.null(reward)){
      if (debug) print("calculating reward")
      cres = calculate_reward(
        formula, 
        id, 
        w = w, 
        tau = tau, 
        data = data,
        owl_method = "MSOWL"
      )
    } else {
      if (debug) print("using precalculated reward")
      cres = reward
    }
    
    # observed data
    dat = cres[[1]] # input data
    complete_data = cres[[2]] # weighted-outcome
    feat = cres[[3]] # features (covariates)
    # print("msowl366:dat")
    # print(dat)
    
    X <- as.matrix(dat[,feat]) # covariates 
    A <- dat$a                 # observed treatment 
    Wt <- complete_data$Wt     # outcome weights 
    
    
    # if (debug) print(Wt)
    
    # compute individualised treatments from rule
    if (class(rule) == "msowl"){
      if (rule$kernel == "linear"){
        form <- as.formula(paste("~", strsplit(as.character(rule$call), "~", fixed = TRUE))[[3]])
        # if (debug) print(form)
        Z <- model.matrix(form, dat)
        f_n <- Z%*%rule$beta_opt
        # if (debug) cat("f_n", f_n)
        
      } else if (rule$kernel == "rbf"){
        rbf <- rbfdot(sigma = rule$sigma)
        # kernel, new_covariates, training_covariates
        K <- kernelMatrix(rbf, as.matrix(dat[,feat]), rule$fit$Z)
        f_n <- rule$fit$beta0 + K%*%rule$fit$alpha1
      }
    } else if(class(rule) == "data.frame"){
      f_n <- rule$A_hat # data.frame [id, A_hat]
    } else{
      # print(paste("owl_method", owl_method))
      print("ERROR ITR NOT OF KNOWN CLASS")
      print(class(rule))
      return(NA)
    }
    
    # Value: compute the empirical expectation (mean)
    # of inverse-probability-of-censoring-weighted outcome 
    # for patients whose observed treatment matches 
    # the prescribed treatment from the itr estimator.
    V_dn <- mean(Wt*(A*f_n>=0))
    # if (debug) cat("V_dn", V_dn)
    
    # Fixed One-Size-Fits-All Rules
    if(fixed.rules){
      f <- 1
      V_1 <- mean(Wt*(A*f>=0))
      f <- -1
      V_m1 <- mean(Wt*(A*f>=0))
      res <- c(V_dn, V_1, V_m1)
      names(res) <- c("V(dn)", "V(1)", "V(-1)")
    } else {
      res <- V_dn
      names(res) <- "V(dn)"
    }
    return(res)
  }
  
  # Compute Value of Treatment Rule on Population Sample
  est <- V_d(data = data)
  
  # Compute Standard Error of Value Estim. via Bootstrap procedure 
  if(SE){
    
    if (debug) {
      pb <- txtProgressBar(
        min = 0,      # Minimum value of the progress bar
        max = nboot,  # Maximum value of the progress bar
        style = 3,    # Progress bar style (also available style = 1 and style = 2)
        width = 50,   # Progress bar width. Defaults to getOption("width")
        char = "="    # Character used to create the bar
      )
    }
    if(!is.character(data$id)){
      data$id <- as.character(data$id)
    }
    clusters <- sort(unique(data[,id]))
    bres <- matrix(NA, nrow=nboot, ncol=(fixed.rules*3 + !fixed.rules*1))
    for(b in 1:nboot){
      index <- sample(1:length(clusters), length(clusters), replace=TRUE)
      cl <- clusters[index]
      bdat <- NULL
      for(j in 1:length(cl)){
        aa <- data[data[,id] %in% cl[j],]
        reps <- table(cl[1:j])[names(table(cl[1:j]))==cl[j]] - 1
        if(reps > 0){
          aa[,id] <- paste(aa[1,id], 1*reps, sep = ".")
        }
        bdat <- rbind(bdat, aa)
      }
      bVal <- try(V_d(data = bdat), silent = TRUE)
      if(class(bVal) != "try-error"){
        bres[b,] <- bVal
      }
      if (debug) {
        setTxtProgressBar(pb, b)
      }
    }
    if (debug) {
      close(pb) # Close the connection
    }

    if(fixed.rules){
      nfail <- sum(is.na(bres[,1]))
      SE <- apply(bres, 2, sd, na.rm = TRUE)

      res <- cbind(est, SE, est - qnorm(0.975)*SE, est + qnorm(0.975)*SE)
      colnames(res) <- c("estimate", "SE", "ll", "ul")
      res <- round(res, 3)

      C <- rbind(c(1, -1, 0),
                 c(1, 0, -1))
      est <- C%*%est
      SE <- sqrt(diag(C%*%var(bres, na.rm = TRUE)%*%t(C)))
      res2 <- cbind(est, est - qnorm(0.975)*SE, est + qnorm(0.975)*SE)
      colnames(res2) <- c("estimate", "ll", "ul")
      rownames(res2) <- c("V(dn) - V(1)", "V(dn) - V(-1)")
      res2 <- round(res2, 3)

      res <- list(res, res2, n.failed = nfail)
    } else {
      nfail <- sum(is.na(bres[,1]))
      SE <- sd(bres, na.rm = TRUE)
      res <- c(est, SE, est - qnorm(0.975)*SE, est + qnorm(0.975)*SE)
      names(res) <- c("estimate", "SE", "ll", "ul")
      res <- list(round(res, 3), n.failed = nfail)
    }
  } else {
    res <- est
  }
  return(res)
}

leave_one_out <- function(i, dat, Xmat, A, Wt, lambda, debug) {
  timer <- Sys.time()
  # print("i")
  # print(i)
  # print("dat")
  # print(dat)
  # print("Xmat")
  # print(Xmat)
  # print("Wt")
  # print(Wt)
  # print("lambda")
  # print(lambda)
  # print("debug")
  # print(debug)
  
  fit1 <- wsvm_solve(
    X = Xmat[dat$id != i,], 
    A = A[dat$id != i], 
    wR = Wt[dat$id != i],
    kernel = 'linear', 
    lambda = lambda
  )
  # if (debug) print(fit1)
  if (unlist(gregexpr('error', class(fit1)))[1] == -1) {
    pred = fit1$beta0 + Xmat[dat$id == i,] %*% fit1$beta
  } else {
    pred = NA
  }
  timer <- Sys.time() - timer
  # if (debug) print(paste("i: ", i, "pred: ", pred, "timer: ", timer))
  return(pred)
}


msowl <- function(
    formula,
    id,
    w = c(0, 1, 0),
    tau = 3,
    data,
    lambda = 1,
    kernel="linear",
    sigma=5,
    jackknife = FALSE,
    jackknife_pg = FALSE,
    trim = NULL,
    owl_method = "MSOWL",
    reward = NULL,
    debug=T
){
  
  if (debug){
    cat("INPUTS: id: ", id, "w: ", w, "tau: ", tau, "dim(data): ", dim(data), "lambda: ", lambda, "kernel:", kernel)
    print(formula)
  }
  
  ### calc_reward ---------
  if (is.null(reward)){
    # if (debug) print("precalc_reward")
    cres = calculate_reward(
      formula, 
      id, 
      w = w, 
      tau = tau, 
      data = data,
      owl_method = owl_method
    )
  } else {
    if (debug) print("using precalculated reward")
    cres = reward
  }
  
  dat = cres[[1]] 
  complete_data = cres[[2]]
  feat = cres[[3]]
  
  X <- as.matrix(dat[,feat])
  A <- dat$a
  Wt <- complete_data$Wt
  
  ### wsvm_solve ------------
  # obtain \hat{f}(x) by minimizing 
  # sum_i=1^n W_i \frac{\phi(A_i f(X_i))}{\pi(A_i;X_i)} + \lambda_n ||f||^2
  fit <- wsvm_solve(
    X=X, 
    A=A, 
    wR=Wt, 
    kernel=kernel, 
    sigma=sigma, 
    lambda=lambda
  )

  # exit if error 
  if(unlist(gregexpr('error', class(fit)))[1] !=-1){
      res <- NA
      class(res) <- "error"
      if (debug){
        print("error!")
        print(fit) 
      }
      return(res)
  }
  
  ### predict ---------
  # assign treatment for the given population from learned itr
  # select subset of patients whose observed treatment matches the ITR assignment
  # compute value of ITR: weighted expectation of survival outcomes.
  if (kernel == "linear"){
    f_n <- fit$beta0 + X%*%fit$beta
    V_dn <- mean(Wt*(A*f_n>=0))
    
    # format output 
    beta_opt <- c(fit$beta0, fit$beta)
    names(beta_opt) <- c("(constant)", feat)
    res <- list(
      beta_opt = beta_opt, 
      Value = V_dn, 
      fit = fit, 
      call = formula,
      kernel=kernel
    )
    
  } else if (kernel=="rbf"){
    # compute treatments for given population from learned itr 
    rbf <- rbfdot(sigma = sigma)
    K <- kernelMatrix(rbf, as.matrix(dat[,feat]), fit$Z)
    f_n <- fit$beta0 + K%*%fit$alpha1
    V_dn <- mean(Wt*(A*f_n>=0))
    
    # format output 
    res <- list(
      Value = V_dn, 
      fit = fit ,
      call = formula,
      kernel=kernel,
      sigma=sigma
    )
  }

  if(kernel == "linear" & isTRUE(jackknife)){
    if (debug) print("jackknife")
    ids <- sort(unique(dat$id))
    f <- rep(NA, times=nrow(dat))
    
    if (jackknife_pg) {
      pb <- txtProgressBar(
        min = 0,      # Minimum value of the progress bar
        max = length(ids),  # Maximum value of the progress bar
        style = 3,    # Progress bar style (also available style = 1 and style = 2)
        width = 50,   # Progress bar width. Defaults to getOption("width")
        char = "="   # Character used to create the bar
      )
    }
    for(i in ids){
      fit1 <- wsvm_solve(
        X=X[dat$id!=i,], 
        A=A[dat$id!=i], 
        wR=Wt[dat$id!=i],
        kernel='linear', 
        lambda=lambda
      )
      if(unlist(gregexpr('error', class(fit1)))[1]==-1){
        f[dat$id==i] <- fit1$beta0 + X[dat$id==i,]%*%fit1$beta
      }
      if (jackknife_pg) setTxtProgressBar(pb, i)
    }
    if (jackknife_pg) close(pb) # Close the connection
    
    ### Sean code added. 
    # Issue with wsvm_solve giving NAs with some jackknife samples (add NA to output)
    vecVal=Wt*(A*f>=0)
    V_jack <- mean(vecVal,na.rm=TRUE)
    SD_jack=sd(vecVal,na.rm=TRUE)
    nMiss_jack=sum(is.na(vecVal))
    res <- list(
      beta_opt = beta_opt, 
      Value = V_dn, 
      Value_jack = V_jack,
      SD_jack=SD_jack,
      nMiss_jack=nMiss_jack, 
      fit = fit, 
      call = formula
    )
  }
  class(res) <- "msowl"
  return(res)
}

# # Select the optimal lambda through leave-one-out cross-validation
select.lambda <- function(
    formula,
    id,
    w,
    tau,
    data,
    lambdas = c(0.01, 0.1, 0.5, 1, 5, 10, 50, 100),
    owl_method = "MSOWL"
){
  V_jack <- rep(NA,length(lambdas))
  SD_jack <- rep(NA,length(lambdas))
  nMiss_jack <- rep(NA,length(lambdas))

  for(i in 1:length(lambdas)){
    print(paste("lambda", lambdas[i], sep = " = "))

    fit <- msowl(
      formula,
      id = id,
      w = w,
      tau = tau,
      data = data,
      lambda = lambdas[i],
      jackknife = TRUE,
      owl_method = owl_method
    )

    if(unlist(gregexpr('error', class(fit)))[1]==-1){
      cat("val", fit$Value_jack, "SD", fit$SD_jack, "nMiss", fit$nMiss_jack)
      V_jack[i] <- fit$Value_jack
      SD_jack[i] <- fit$SD_jack
      nMiss_jack <- fit$nMiss_jack
    } else{
      print("error ... :( ")
    }
  }
  #Should pick lambda better with standard deviation as well
  lambda <- min(lambdas[V_jack==max(V_jack, na.rm = TRUE)], na.rm = TRUE)
  res <- list(
    lambda = lambda,
    details = cbind(lambdas, V_jack,SD_jack,nMiss_jack)
  )
  return(res)
}
# select.lambda <- function(
#   formula, 
#   id, 
#   w, 
#   tau, 
#   data, 
#   lambdas = c(0.01, 0.1, 0.5, 1, 5, 10, 50, 100),
#   owl_method = "MSOWL",
#   debug = F,
#   jackknife_pg = F
# ){
#   V_jack <- rep(NA,length(lambdas))
#   SD_jack <- rep(NA,length(lambdas))
#   nMiss_jack <- rep(NA,length(lambdas))
#   for(i in 1:length(lambdas)){
#     if (debug) print(paste("lambda", lambdas[i], sep = " = "))
#     fit <- msowl(
#       formula, 
#       id = id, 
#       w = w, 
#       tau = tau, 
#       data = data,
#       lambda = lambdas[i] * length(unique((data$id)))^(-1/2), 
#       jackknife = TRUE,
#       jackknife_pg = F,
#       owl_method = owl_method,
#       debug = debug
#     )
#     # print(fit)
#     if(unlist(gregexpr('error', class(fit)))[1]==-1){
#       V_jack[i] <- fit$Value_jack
#       SD_jack[i] <- fit$SD_jack
#       nMiss_jack <- fit$nMiss_jack
#     }
#   }
#   #Should pick lambda better with standard deviation as well
#   lambda <- min(lambdas[V_jack==max(V_jack, na.rm = TRUE)], na.rm = TRUE)
#   res <- list(lambda = lambda, details = cbind(lambdas, V_jack,SD_jack,nMiss_jack))
#   return(res)
# }

#### Examples ####
# select.lambda(ms(t1 = t1, t2 = t2, s1 = s1, s2 = s2, a = A) ~ Z1 + Z2,
#                 id = "id", w = c(0, 1, 0), tau = 3, data = dat,
#                 lambdas = c(0.01, 0.1, 0.5, 1, 5, 10, 50, 100)
# )
# 
# fit <- msowl(ms(t1 = t1, t2 = t2, s1 = s1, s2 = s2, a = A) ~ Z1 + Z2,
#             id = "id", w = c(0, 1, 0), tau = 3, data = dat,
#             lambda = 1
# )
# Value <- msowl.val(
#   ms(t1 = t1, t2 = t2, s1 = s1, s2 = s2, a = A) ~ Z1 + Z2,
#   id = "id",
#   w = c(0, 1, 0),
#   tau = 3,
#   rule = fit,
#   fixed.rules = TRUE,
#   data = dat
# )




cross.validation <- function(
    formula,
    data,
    id = "id",
    w = c(1,0),
    tau = 10,
    lambda = 1,
    k = 5,
    debug = F,
    SE=F,
    owl_method = "MSOWL",
    B = 10
  ){
  timer = Sys.time()
  V     <- rep(NA,B)
  SD    <- rep(NA,B)
  if (debug){
    pb <- txtProgressBar(
      min = 0,      # Minimum value of the progress bar
      max = B*k,  # Maximum value of the progress bar
      style = 3,    # Progress bar style (also available style = 1 and style = 2)
      width = 50,   # Progress bar width. Defaults to getOption("width")
      char = "="    # Character used to create the bar
    )
  }
  # B REPETITIONS
  for (i in 1:B) {
    if (debug) print(paste("b = ", i))
    # RANDOM K SUBSETS 
    folds = sample.int(
      n=k,
      size=nrow(data),
      replace=T
    )
    # per round results
    cox.V.k     <- rep(NA,k)
    V.k     <- rep(NA,k)
    SD.k    <- rep(NA,k)
    for (j in 1:k){
      if (debug) {print(paste0("fold: ", j, "/", k) )}
      train_set = as.data.frame(data[which(folds != j), ])
      test_set  = as.data.frame(data[which(folds == j), ])
      if (owl_method %in% c("MSOWL", "ICO")) {
        fit <- msowl(
          formula,
          id = "id",
          w = w,
          tau = tau,
          data = train_set,
          lambda =  nrow(train_set)^(-1/2)*lambda,
          jackknife = FALSE,
          owl_method = owl_method
        )
        if(unlist(gregexpr('error', class(fit)))[1]==-1){
          val <- msowl.val(
            formula=formula,
            id="id",
            w=w,
            tau=tau,
            rule = fit,
            data = test_set,
            fixed.rules=T,
            SE = SE
          )
          if (debug) print(val)
          # aggregate results
          if (SE){
            V.k[j] <- val[[1]][1,1] # val
            SD.k[j] <- val[[1]][1,2] # SE
          } else{
            V.k[j] <- val[[1]][1] # val
          }
          
        } else {
          if (debug) print("error ... :( ")
          V.k[i] <- NA
          SD.k[i] <- NA
        }
        
      }
      # COX -----------
      print("COX")
      # FIT MODEL
      pfs.cox.formula = as.formula("Surv(time=pfs, pfs.event.status, type='right') ~ age + factor(sex) + factor(ecog) + factor(kras) + factor(a) + factor(a):age + factor(a):factor(sex) + factor(a):factor(ecog) + factor(a):factor(kras)")
      cox.train <- coxph(
        formula = pfs.cox.formula, 
        data = train_set
      )
      if (debug) summary(cox.train)
      # test.ph <- cox.zph(cox.train)
      # test.ph
      # ggcoxzph(test.ph)
      # EVAL TEST SET
      # cox.test.predictions <- coxph_itr_predict(
      #   cox_fit = cox.train,
      #   data = test_set
      # )
      data.copy = copy(test_set)
      data.copy$a = 1
      res = survfit(cox.train, newdata= data.copy)
      pred_pos = summary(res)$table[,"median"]
      # -1
      data.copy$a = -1  
      res = survfit(cox.train, newdata= data.copy)
      pred_neg = summary(res)$table[,"median"]
      if (any(is.na(pred_pos)) | any(is.na(pred_pos))) {
        print("NA in PREDICTION! ")
        print(head(data))
        print("NA in PREDICTION! ")
      }
      better_treatment = ( which.pmax(pred_neg, pred_pos) - 1 ) * 2 - 1
      cox.test.predictions = data.frame(test_set$id, A_hat=better_treatment)
      
      val <- msowl.val(
        formula=formula,
        id="id",
        w = c(1,0),
        tau=tau,
        rule = cox.test.predictions,
        fixed.rules = T,
        data = test_set
      )
      if (debug) print(val)
      # TODO aggregrate result
      cox.V.k[j] = val[1]
    
      if (debug) setTxtProgressBar(pb, (i-1)*k+j-1)
    }
    if (debug) close(pb) # Close the connection
    cox.V <- mean(cox.V.k, na.rm=T)
    V[i] <- mean(V.k,na.rm=T)
    SD[i] <- mean(SD.k,na.rm=T)
  }
  
  res = list(
    V = V,
    SD = SD,
    cox.V = cox.V,
    duration =  Sys.time() - timer
  )
  return(res)
}


# Utility Functions  -----------------------------------------------------------
## K-Fold CV SVM Tuning --------------------------------------------------------
select.lambda.kfold <- function(
    formula, 
    id, 
    w, 
    tau, 
    data, 
    lambdas = c(0.01, 0.1, 0.5, 1, 5, 10, 50, 100),
    k = 5,
    debug = F,
    SE=T,
    owl_method = "MSOWL",
    round_vals = NULL
){
  # print(owl_method)
  timer = Sys.time()
  V     <- rep(NA,length(lambdas))
  SD    <- rep(NA,length(lambdas))
  
  if (debug){
    pb <- txtProgressBar(
      min = 0,      # Minimum value of the progress bar
      max = length(lambdas)*k,  # Maximum value of the progress bar
      style = 3,    # Progress bar style (also available style = 1 and style = 2)
      width = 50,   # Progress bar width. Defaults to getOption("width")
      char = "="    # Character used to create the bar
    )
  }
  for(i in 1:length(lambdas)){
    if (debug) print(paste0(i,") lambda = ", lambdas[i]))
    folds = sample.int(
      n=k,
      size=nrow(data), 
      replace=T
    )
    V.k     <- rep(NA,k)
    SD.k    <- rep(NA,k)
    for (j in 1:k){
      if (debug) {
        print(paste("fold: ", j))
      }
      fit <- msowl(
        formula,
        id = "id", 
        w = w, 
        tau = tau, 
        data = as.data.frame(data[which(folds != j), ]),
        lambda =  nrow(data[which(folds != j), ])^(-1/2)*lambdas[i], 
        jackknife = FALSE,
        owl_method = owl_method
      )
      if(unlist(gregexpr('error', class(fit)))[1]==-1){
        val <- msowl.val(
          formula=formula,
          id="id",
          w=w,
          tau=tau,
          rule = fit,
          data = data[which(folds ==j),],
          # lambda=lambdas[i],
          fixed.rules=T,
          SE = SE
        )
        if (debug) print(val)
        if (SE){
          V.k[j] <- val[[1]][1,1] # val
          SD.k[j] <- val[[1]][1,2] # SE 
        } else{
          V.k[j] <- val[[1]][1] # val
        }
      } else{
        if (debug) {
          print("error ... :( ")
        }
        V.k[i] <- NA
        SD.k[i] <- NA
      }
      if (debug){
        setTxtProgressBar(pb, (i-1)*k+j-1) 
      }
    }
    if (debug){
      close(pb) # Close the connection  
    }
    V[i] <- mean(V.k,na.rm=T)
    SD[i] <- mean(SD.k,na.rm=T)
    
    if (debug) {
      print("value per k-folds for this lambda: ")
      print("V.k")
      print(V.k)
      if (SE) {
        print("SD")
        print(SD.k)
      }
      print("lambda avg values until now: ")
      print(V)
      if (SE) {
        print("SD")
        print(SD)
      }
    }
  }
  
  # if (round_vals == "1se"){
  # get lambda with minimum value that is within 1se of the max val 
  # details = setDT(as.data.frame(pfs.lambda$details))
  # max_val = details$V == max(details$V, na.rm=T)
  # lambda.max <- details[max_val,]
  # lambda.1se = details[V + 2*SD >= lambda.max$V, min(lambdas)]
  
  # 1.  pure maximum (with min variance tie-breaker)
  lambda.max <- min(lambdas[V==max(V, na.rm = TRUE)], na.rm = TRUE)
  
  # 2. Rounded lambda max  
  # V_rounded = sapply(FUN=round, X=V, round_vals)
  # lambda.max.rounded <- min(lambdas[V_rounded==max(V_rounded, na.rm = TRUE)], na.rm = TRUE)
  
  res <- list(
    lambda.max = lambda.max, 
    # lambda.max.rounded = lambda.max.rounded, 
    lambdas = lambdas,
    # details = cbind(lambdas, V, SD),
    details = cbind(lambdas, V, SD),
    duration = Sys.time() - timer,
    owl_method = owl_method,
    k = k,
    SE = SE
  )
  
  # 3. lambda.1se: min lambda whose value is within 1se of max val
  details = setDT(as.data.frame(res$details))
  max_val = details$V == max(details$V, na.rm=T)
  lambda.max <- details[max_val,]
  lambda.1se = details[V + SD >= lambda.max$V, min(lambdas)]
  res$lambda.1se = lambda.1se
  
  if (debug) {
    print("FINAL")
    print(res)
  }
  return(res)
}


precalculate_reward = function(
    data,
    # n,
    formula,
    tau,
    w,
    id,
    n.batch,
    owl_method="MSOWL",
    debug=F
) {
  if (debug) print("pre-calculate reward in small batches")
  n_groups = nrow(data) / n.batch
  data.batches = split(
    as.data.frame(data),
    sample.int(
      n=n_groups,
      size=nrow(data),
      replace=TRUE
    )
  )
  cres = pbapply::pbsapply(
    # cres = sapply(
    X=data.batches,
    FUN=calculate_reward,
    formula = formula,
    id = id,
    w = w,
    tau = tau,
    owl_method=owl_method
  )
  res = list(
    dat = bind_rows(cres[1,]),
    complete_data = bind_rows(cres[2,]),
    feat = cres[3,][[1]]
  )
  res$dat = res$dat[order(res$dat$id), ]
  res$complete_data = res$complete_data[order(res$complete_data$id), ]
  return(res)
}

# Divide & Conquer Value Function -------------------------------------------
msowl.itr.divide.conquer.value = function(
    data,
    rule,
    formula,
    id="id",
    w = c(1,0),
    tau=tau,
    fixed.rules=T,
    # n.test,
    n.batch,
    jackknife = F
){
  batches = split(
    as.data.frame(data),
    sample.int(
      n=nrow(data)/n.batch,
      size=nrow(data),
      replace=TRUE
    )
  )
  res = pbapply::pbsapply(
    batches,
    msowl.val,
    formula=formula,
    id="id",
    w = w,
    tau=tau,
    rule=rule,
    fixed.rules=fixed.rules
  )
  if (fixed.rules){
    res = c(
      "V(dn)" = mean(res[1,]),
      "V(-1)" = mean(res[2,]),
      "V(1)" = mean(res[3,])
    )
  } else {
    mean(res["V(dn)"])
  }
  return(res)
}

# DIVIDE & CONQUER for COX Regression Model ------------------------------------
itr.cox.divide.conquer.value = function(
    data,
    cox_fit,
    formula,
    id="id",
    w = c(1,0),
    tau=tau
){
  cox.test.predictions <- coxph_itr_predict(
    cox_fit = cox_fit,
    data = data,
    tau=tau
  )
  
  cval = msowl.val(
    data=data,
    formula=formula,
    id="id",
    w = c(1,0),
    tau=tau,
    rule=cox.test.predictions,
    fixed.rules=F
  )
}

# cox.test.val = pbapply::pbsapply(
#   X=data.test.batches,
#   FUN=itr.cox.divide.conquer.value,
#   cox_fit = cox.train,
#   formula=my_formula,
#   id="id",
#   w = c(1,0),
#   tau=tau
# )

