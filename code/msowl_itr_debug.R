########
# Estimation of an optimal individualized treatment rule
library(survival)
library(ggplot2)
library(ggfortify)
library(ranger)
library(dplyr)

# 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#########################################
# DATA: SPECTRUM OS (Overall Survival)
#########################################
spectrum_os = read.csv("~/aueb/dtr/data/spectrum_transformed_data/spectrum_os.csv", sep = "\t")
spectrum_os$censor_status = spectrum_os$s1 == spectrum_os$s2 
spectrum_os = mutate(spectrum_os, A = factor(A,labels=c("C", "C+P")))
# kaplan meier curve
km = with(spectrum_os, Surv(t2, censor_status))
km_fit = survfit(
  Surv(t2,censor_status) ~ A, 
  data=spectrum_os
)
summary(km_fit)
autoplot(km_fit)
# ggsave("spectrum_os_km_curve.pdf")
ggsave("~/aueb/dtr/documents/thesis_paper/diagrams/spectrum_os_km_curve.pdf", width=10, height=5)



#########################################
# DATA: SPECTRUM MS (Multistate) 
#########################################
spectrum_ms = read.csv("~/aueb/dtr/data/spectrum_transformed_data/spectrum_multistate.csv", sep = "\t")

# censored status: if record has same start and end state
spectrum_ms$censor_status = spectrum_ms$s1 == spectrum_ms$s2

# pretty label for plot
spectrum_ms = mutate(spectrum_ms, A = factor(A,labels=c("C", "C+P")))

# state 1 
sms_s1 = spectrum_ms[spectrum_ms$s1 == 1, ]
# kaplan meier curve
km = with(sms_s1, Surv(t2-t1, s1 == s2))
km_fit = survfit(
  Surv(t2-t1, s1 == s2) ~ A, 
  data=sms_s1
)
summary(km_fit)
autoplot(
  km_fit, 
  xlim=c(0,540), 
  ylim=c(0.65,1), 
  main="Initial State"
)

ggsave("~/aueb/dtr/documents/thesis_paper/diagrams/spectrum_ms_s1_km_curve.pdf", width=10, height=5)


# state 2 
sms_s2 = spectrum_ms[spectrum_ms$s1 == 2, ]
# kaplan meier curve
km = with(sms_s2, Surv(t2-t1, s1 == s2))
km_fit = survfit(
  Surv(t2-t1, s1 == s2) ~ A, 
  data=sms_s2
)
summary(km_fit)
autoplot(
  km_fit, 
  xlim=c(0,540), 
  ylim=c(0.65,1), 
  main="Progressive Disease State"
)

# ggsave("spectrum_ms_s2_km_curve.pdf")
ggsave("~/aueb/dtr/documents/thesis_paper/diagrams/spectrum_ms_s2_km_curve.pdf", width=10, height=5)


####################################
#
###################################

sim1 = read.csv("/Users/andreehrlich/aueb/dtr/data/simulated_data_n=500_scenario=1_cr=-1.6.csv",sep="\t")
sim1

#########################################
# 
#########################################

# itr <- function(
  
#### spectrum OS 
# data = spectrum_os
# feat = c("Z1", "Z2", "Z3", "Z4")
# w=c(1) # c(0, 1, 0) 
# tau=540 
# kernel='linear' 
# lambda=1
# sigma=1 
# SE=TRUE

### simulated multi_state
data = sim1
feat = c("Z1", "Z2") #, "Z3", "Z4")
w=c(1, 1, 0) 
tau=3 
kernel='linear' 
lambda=1
# sigma=1 
SE=TRUE
  
  # ){

    if(!(kernel %in% c('linear', 'rbf'))){
      stop("Only 'linear' and 'rbf' (radial basis function or Gaussian) kernels are supported")
    }
    
    # Multi-State formating
    # state space
    S <- with(data, sort(unique(c(s1, s2))))

    ## absorbing state subspace
    T <- S[!(S %in% sort(unique(data$s1)))]

    if(length(S)!=length(w)){
      stop("length(w) not equal to the total number of states")
    }

    if(min(w)<0 | max(w)>1 | sum(w)>=length(w)){
      stop("min(w)<0 or max(w)>1 or sum(w)>=length(w)")
    }

    # Survival Time Formating  
    if(max(data$t1>data$t2)==1){
      stop("t1 > t2")
    }
    
    # data[which(data$t1>data$t2),]

    # Calculate Reward 
    # Estimate censoring distribution 

    # censored status
    c.delta <- 1*(data$s1==data$s2 & data$t2<tau)
    
    # utilise coxph model just for the built-in 
    # non-parametric Nelson-Aalen estimator
    fit <- coxph(
      Surv(
        time=t1, 
        time2=t2, 
        event=c.delta,
        type="counting")
      ~ 1, 
      data=data
    )
    # summary(fit)
    # predicted survival curve for cox model 
    # instantaneous rate of death computed at points in time at which death is observed 
    Haz <- basehaz(fit, centered=FALSE) 
    str(Haz$time)
    Haz
    # get uncensored event times 
    uncensored = data[data$s1!=data$s2,"t2"] # uncensored 
    str(uncensored)
    uncensored
    
    # tt: all of the merge predicted and actual event times
    tt <- sort(unique(c(Haz$time, uncensored)))
    tt <- c(0,tt[tt<=tau]) # non-censored times only ! 
    
    # zeta_t function 
    # 
    event <- matrix(
      sapply(tt, zeta_t, data=data, w=w, hazard=Haz, S=S, T=T), 
      byrow=TRUE, 
      nrow=length(tt)
    )
    
    # integrate inverse-probability-of-censoring weighted (IPCW) outcomes over time (in discrete states)
    # time between events 
    dm <- diff(c(tt, tau), lag=1)
    xi <- colSums(event*dm).     
    xi <- matrix(tapply(X=xi, INDEX = data$id, FUN=sum))
    
    # then do inverse-propensity-of-treatment (IPW )
    # propensity score 
    dat <- data[data$t1==0,]
    dat <- dat[order(dat$id),]
    pi_n <- mean(dat$A==1)
    
    # weighting 
    Wt <- xi/(pi_n*(dat$A==1) + (1 - pi_n)*(dat$A==-1))
    
    # patient covariates
    X <- as.matrix(dat[,feat])
    
    # treatment 
    A <- dat$A
    
    
    #  weighted SVM 
    # Note: Currently hardcoded to choose either linear or rbf kernel
    if(kernel=='linear'){
      fit <- wsvm_solve(X=X, A=A, wR=Wt, kernel='linear', lambda=lambda)
      
      if(!is.na(fit$beta0)){
        beta_opt <- c(fit$beta0, fit$beta)
        V_opt <- V_d(data=data, w=w, tau=tau, dec.fun=fit, feat=feat, SE=SE)
        res <- list(beta_opt=beta_opt, V_opt=V_opt$V_n, se_V_opt=V_opt$se, fit=fit)
      } else {
        res <- NA
        class(res) <- "error"
      }
    } else {
      fit <- wsvm_solve(X=X, A=A, wR=Wt, kernel='rbf', sigma=sigma, lambda=lambda)
      if(!is.na(fit$beta0)){
        V_opt <- V_d(data=data, w=w, tau=tau, dec.fun=fit, feat=feat, SE=FALSE)
        res <- list(V_opt=V_opt, fit=fit)
      } else {
        res <- NA
        class(res) <- "error"
      }
    }
    return(res)
  }
    