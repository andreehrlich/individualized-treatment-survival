# Survival Regression Models 
# - COX 
# - AFT 

# For eeach, 
# - Generate Survival Times 
# - ITR Outcome Regression 

# Simulate survival times from cox regression 
# Simulate Data from a Cox Regression  DONE
# (P = 1,2,3 var + interaction terms).    
# Fit cox model (full, main fx only, intx-fx only). DONE
# Fit 100x and take mean/std of parameter estimates. DONE 
# Compare results given different baseline hazard functions  TODO 


# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# getwd()
library(data.table)
library(survival)
library(survminer)
library(radiant.data)

# survival generation packages
library(flexsurv)
library(eha)


source("./msowl_modified.R")


# AFT --------
# Accelerated Failure Time Regression: Log-Normal
sim.surv.regr <- function(n, H_inv, psi, model)  {
  U <- runif(n)
  if (model == "aft") {
    times = exp(-log(1-U)) / exp(psi) # log-normal baseline
  } else if (model == "cox") {
    times = exp(qnorm(U)) * exp(-psi) # mu = 0, sigma = 1
  }
  return(times)
}

# BAS HAZ FUNC -------
# TODO try with bas.haz distributions: (1) exponential & (2) weibull
# kosorok2015 paper: bas_haz = function(t) 2*t  => cum_bas_haz = function(t) t^2 
H.inv.kos =function(t) t^(1/2)
# COX ITR --------
coxph_itr_predict <- function(cox_fit, data){
  data.copy = copy(data)
  data.copy$A = 1
  res = survfit(cox_fit,newdata= data.copy)
  pred_pos = summary(res)$table[,"rmean"]
  # -1
  data.copy$A = -1  
  res = survfit(cox_fit,newdata= data.copy)
  pred_neg = summary(res)$table[,"rmean"]
  if (any(is.na(pred_pos)) | any(is.na(pred_pos))) {
    print("NA in PREDICTION! ")
    print(head(data))
    print("NA in PREDICTION! ")
  }
  better_treatment = ( which.pmax(pred_neg, pred_pos) - 1 ) * 2 - 1
  itr = data.frame(data$id, A_hat=better_treatment)
  return(itr)
}

# Data Generation ----------------------------------------------------------------------
# ZhaoKosorok2015 OWL Methods for Right-Censored Data Simulation Replication
# scenario.vec = 1:10
# tau.vec.OLD = c(1.5, 2, 2, 2.5, 3, 3, 3, 3, 2, 2)
tau.vec = c(1.5, 2, 2, 2.5)
# tau.vec = tau.vec.KOS + 2
sim.kos = function(
    n = 500,
    scenario = 1,
    tau = tau.vec[1],
    censor.rate = -1
){
  A = rbinom(n=n, 1, prob = 0.5) * 2 - 1 
  z1 = runif(n,0,1)
  z2 = runif(n,0,1)
  z3 = runif(n,0,1)
  # CENSORING
  if (censor.rate == FALSE) {
    time.vec.censoring = rep(10000,n) 
  } else {
    time.vec.censoring = rexp(n,rate= exp(censor.rate)) 
  }
  # SURVIVAL
  if (scenario == 1){
    time.vec.survival = sim.surv.regr(
      n = n, 
      H_inv = H.inv.kos, 
      psi = 0.6*z1 - 0.8*z2 + (0.6 - 0.4*z1 - 0.2*z2 - 0.4*z3)*A,
      model = "cox"
    )
  } else if (scenario == 2){
    time.vec.survival = sim.surv.regr(
      n=n, 
      H_inv=H.inv.kos, 
      psi = -1.5*z1 + 0.5*z2 + (1 - 0.7*z1 - 1.2*z2)*A,
      model = "cox"
    )
  } else if (scenario == 3){
    time.vec.survival = sim.surv.regr(
      n = n, 
      H_inv = H.inv.kos, 
      psi = -0.5 - 0.8*z1 + 0.7*z2 + 0.2*z3 +(0.6-0.4*z1 -0.1*z2-0.4*z3)*A,
      model = "aft"
    )
  } else if (scenario == 4){
    time.vec.survival = sim.surv.regr(
      n = n, 
      H_inv = H.inv.kos, 
      psi = -0.2 - 0.5*z1 + 0.5*z2 + 0.3*z3 +(0.5-0.1*z1 -0.6*z2+0.1*z3)*A, 
      model = "aft"
    )
  } else if (scenario %in% c(5,6,7,8)) {
    # bakoyannis scenarios: transition intensity for treatment effect function 
    f_star = function(z1, z2, scenario) {
      switch(as.character(scenario),
             `1` = z1 + z2,
             `2` = 2 * (z1 - z2),
             `3` = 1 + z2 - exp(-z1),
             `4` = 2 * log(2 - z1 - z2) - 1.4)
    }
    # randomly sample the theta parameter (mean transition intensity) 
    # and then sample from exponential ? 
    transition_intensity = rexp(n, -0.5 * z1 + 0.5 * z2 + A * f_star(z1, z2, scenario - 4))
    survival_times = rexp(n,transition_intensity) 
  } else if (scenario == 9){
    time.vec.survival = sim.surv.regr(
      n=n, 
      H_inv=H.inv.kos, 
      psi = -1.5*z1 + 0.5*z2 + (1 - 0.7*exp(z1) - 1.2*z2^2 )*A,
      model="cox"
    )
  } else if (scenario == 10){
    time.vec.survival = sim.surv.regr(
      n = n, 
      H_inv = H.inv.kos, 
      psi = -0.2 - 0.5*z1 + 0.5*z2 +(0.5-0.1*z1 -0.6*z2)*A,
      model = "aft"
    )
  }
  # OBSERVED INFORMATION 
  # if (censor.rate == FALSE) {
  #   tau = 10000000000
  # }
  
  # restricted survival time by end-of-study
  time.vec.survival = pmin(time.vec.survival, tau)
  # censoring
  time.vec.observed = pmin(time.vec.survival, time.vec.censoring)
  status = as.numeric(which.pmin(time.vec.survival, time.vec.censoring) == 1)
  # end state 
  state.terminal = status + 1   # state 1 = alive; state 2 = dead
  records = data.frame(
    id = 1:n,
    z1 = z1,
    z2 = z2,
    z3 = z3,
    A  = A,
    s1 = 1,
    s2 = state.terminal,
    t1 = 0,
    t2 = time.vec.observed,
    status = status
  )
  records
  return(records)
}



################
# DATA GENERATION TESTING -------
test_dgf <- function(
  n = 100,
  censor = -1,
  scenario = 2
) {
  tau = tau.vec[scenario] 
  df = sim.kos(
    n=n,
    scenario=scenario,
    tau=tau,
    censor.rate = censor
  )
  df = setDT(df)
  # df[t2 > tau,]
  unique(df[status == 0, s2])
  ggplot(df, aes(x=t2, color=factor(status))) + geom_histogram(bins = n)
  cat(mean(!df$status))

  fit.aft <- aftreg(
    formula = Surv(t2, status) ~  z1 + z2 + z3 + A + A:z1 + A:z2 + A:z3, 
    data=df, 
    dist="lognormal"
  )
  fit.aft
  summary(fit.aft)
  # fit models
  # fit.aft = lm(
  #   formula = log(t2) ~ z1 + z2 + z3 + A + A:z1 + A:z2 + A:z3 , 
  #   data= df
  # )
  fit.cox <- coxph(
    formula= Surv(time=t2, status) ~ z1 + z2 + z3 + A + A:z1 + A:z2 + A:z3,
    data=df
  )
  fit.cox
  summary(fit.cox)
  
  lambda.tuning = select.lambda(
    ms(t1 = t1, t2 = t2, s1 = s1, s2 = s2, a = A) ~ z1 + z2 + z3, 
    id = "id", 
    w = c(1,0), 
    tau=tau, 
    data = as.data.frame(df), 
    lambdas = lambda.vec,
    owl_method = "MSOWL",
    debug = T
  )
  
  lala = as.data.frame(lambda.tuning$details)
  ggplot(lala, aes(x=log(lambdas), y = V_jack)) + 
    geom_point() + 
    geom_point(aes(x=log(lambdas), y=SD_jack, col="red"))
  
  lambda.max_val <- min(lala$lambdas[lala$V_jack==max(lala$V_jack, na.rm = TRUE)], na.rm = TRUE)
  if (!is.null(round_vals)) {
    # V_rounded = sapply(FUN=round, X=V_jack, 1)
    lambda <- min(lambdas[V_rounded==max(V_rounded, na.rm = TRUE)], na.rm = TRUE)
  }
  
  msowl_train <- msowl(
    formula = ms(t1 = t1, t2 = t2, s1 = s1, s2 = s2, a = A) ~ z1 + z2 + z3,
    id = "id",
    w = c(1,0),
    tau = tau,
    data = as.data.frame(df),
    lambda = nrow(df)^(-1/2)*lalalaksdjfalskdfjasldfk, 
    jackknife = T,
    trim = NULL,
    reward = NULL,
    debug=T
  )
  msowl_train$Value
  
  msowl_val = msowl.val(rule = msowl_train, fixed.rules = T, SE=T)
  
  # itr.cox.val = pbapply::pbsapply(
  #   X=batches,
  #   FUN=itr.cox.divide.conquer.value,
  #   cox_fit=cox.train,
  #   formula=formula,
  #   id="id",
  #   w = c(1,0),
  #   tau=tau
  # )
}



# COX ITR SIM -------
cox.itr.sim <- function(
    scenario = 1,
    n.test = 10000,
    n.train = 400,
    tau = 2.5 # tau.vec[scenario]
){
  formula = ms(t1 = t1, t2 = t2, s1 = s1, s2 = s2, a = A) ~ z1 + z2 + z3
  # TRAIN
  data.train = sim.kos(
    n=n.train,
    scenario=1,
    tau=tau.vec[1],
    censor.rate=F
  )
  # TEST 
  data.test = sim.kos(
    n=n.test,
    scenario=scenario,
    tau=tau.vec[scenario],
    censor.rate=F
  )
  data.test.complete = FALSE
  batches = split(
    as.data.frame(data.test),
    sample.int(
      n=n.test/n.batch,
      size=nrow(data.test),
      replace=TRUE
    )
  )
  
  ## PRECALC REWARD --------
  # data.test.complete = precalculate_reward(
  #   data=data.test,
  #   n=10000,
  #   formula=my_formula,
  #   tau=tau.vec[1],
  #   n.batch = 200
  # )
  
  timer = Sys.time()
  print(timer)
  cres = pbapply::pbsapply(
    X=batches,
    FUN=calculate_reward,
    formula = formula,
    id = "id",
    w = c(1,0),
    tau = tau,
    owl_method="MSOWL"
  )
  timer = Sys.time() - timer
  print(timer)
  # print(cres)
  res = list(
    dat = bind_rows(cres[1,]),
    complete_data = bind_rows(cres[2,]),
    feat = cres[3,][[1]]
  )
  res$complete_data = res$complete_data[order(res$complete_data$id),]
  
  
  ## COX -------
  cox_fit <- coxph(
    formula=as.formula("Surv(time=t1, time2=t2, status) ~ z1 + z2 + z3 + factor(A) + factor(A):z1 + factor(A):z2 + factor(A):z3"), 
    data=data.train
  )
  summary(cox_fit)
  # PREDICT ON TEST DATA 
  cox.test.predictions = coxph_itr_predict(
    cox_fit = cox_fit,
    data = data.test,
    tau.vec[1]
  )
  itr.cox.val <- msowl.val(
    formula = formula,
    id = "id",
    w = c(1,0),
    tau = tau,
    rule = cox.test.predictions,
    fixed.rules = FALSE,
    data = as.data.frame(data.test),
    reward = data.test.complete
  )
}

# SIM AFT + COX ----------
## SCEN COEFS  --------
scen.coef = list(
  # SCEN 1 
  setDT(data.frame(
    "variable"= c( "z1", "z2", "z3", "A", "z1:A", "z2:A", "z3:A" ),
    "value" =   c( 0.6,  -0.8,    0, 0.6,   -0.4,   -0.2,  -0.4  )
  )),
  # SCEN 2
  setDT(data.frame(
    "variable"= c( "z1", "z2", "A", "z1:A", "z2:A"),
    "value" =   c( -1.5,  0.5,  1,   -0.7,   -1.2 )
  )),
  # scen 3
  setDT(data.frame(
    "variable"= c( "(Intercept)", "z1", "z2", "z3", "A", "z1:A", "z2:A", "z3:A" ),
    "value" =   c(          -0.5, -0.8,  0.7,  0.2, 0.6,   -0.4,   -0.1,   -0.4 )
  )),
  # SCEN 4
  setDT(data.frame(
    "variable"= c( "(Intercept)", "z1", "z2", "z3", "A", "z1:A", "z2:A", "z3:A" ),
    "value" =   c(          -0.2, -0.5,  0.5,  0.3, 0.5,   -0.1,  -0.6,   0.1 )
  ))
)

## SIM FUNC --------
sim.surv.cox.aft <- function(
    n = 100,
    censor.vec = c(F, -1),
    B = 1000,
    scenario = 3
){
  for (scenario in c(1:4)) {
    for (censor in censor.vec) {
      for (b in 1:B) {
        
        df = sim.kos(
          n=n,
          scenario=scenario,
          tau=tau.vec[scenario],
          censor.rate = censor
        )
        # n = 100
        # censor = F
        # scenario = 1
        # tau = tau.vec[scenario]
        # 
        # fit.aft <- aftreg(
        #   formula = Surv(t2, status) ~  z1 + z2 + z3 + A + A:z1 + A:z2 + A:z3, 
        #   data=df, 
        #   dist="lognormal"
        # )
        
        # fit models
        fit.aft = lm(
          formula = log(time) ~ z1 + z2 + z3 + A + A:z1 + A:z2 + A:z3 , 
          data= df
        )
        fit.cox <- coxph(
          formula= Surv(time=time, status) ~ z1 + z2 + z3 + A + A:z1 + A:z2 + A:z3,
          data=df
        )
        # save res
        if (b == 1) {
          res.aft = matrix(NA, nrow=B, ncol=length(fit.aft$coefficients)) 
          res.cox = matrix(NA,nrow=B, ncol=length(fit.cox$coefficients))
        }
        res.aft[b, ] = fit.aft$coefficients # AFT
        res.cox[b, ] = fit.cox$coefficients # COX
        if (b == B){
          # AFT
          res.aft = setDT(data.frame(res.aft))
          colnames(res.aft) = names(fit.aft$coefficients)
          # COX
          res.cox = setDT(data.frame(res.cox))
          colnames(res.cox) = names(fit.cox$coefficients)
        }
      } # B REPETITIONS
      # PLOT FOR CURRENT ROUND! 
      res.aft$model = "AFT"
      res.cox$model = "COX"
      # GGPLOT HISTOGRAMS
      # AFT
      dt = copy(res.aft)
      mvars = c("(Intercept)", "z1","z2","z3","A","z1:A","z2:A","z3:A")
      # dt[, model:= "AFT"]
      dt.melt = melt(dt, measure.vars = mvars,id.vars = c("model"))
      dt.melt$variable = factor(dt.melt$variable,mvars)
      # COX
      dt2 = copy(res.cox)
      cnames = copy(colnames(dt2))
      # dt2$model = "COX"
      cnames = c("A", "z1", "z2", "z3", "z1:A", "z2:A", "z3:A")
      print(cnames)
      dt2.melt = melt(dt2,measure.vars = cnames, id.vars = c("model"))
      # JOIN TO ONE DT 
      dt3 = rbind(dt.melt,dt2.melt)
      model.colors <- c(
        COX = "orange",
        AFT = "blue"
      )
      # PLOT
      ggplot(
        dt3, 
        aes(x=value, fill=model) # factor(model, color=color))
      ) + 
        scale_fill_manual(values=model.colors) + 
        geom_histogram(
          binwidth = 0.01,
          alpha = 0.75
        ) + 
        facet_wrap(
          ~ variable, ncol=4
        ) +
        geom_vline(
          aes(xintercept=value, color="red"), 
          scen.coef[scenario][[1]]
        ) +
        theme(
          legend.position = NULL # "bottom"
        )
      fname = paste0("./documents/thesis_paper/figures/sim_surv_hist_scen",scenario,".pdf")
      print(paste0("saving plot to disk: ", fname))
      ggsave(fname,width=8,height = 6)
      
      # Put all Statistics around Parameter Estimates in Table 
      dt.table = setDT(dplyr::bind_rows(res.aft, res.cox))
      dt.table = dt.table[, 
        sapply(.SD, function(x) list(paste0(round(mean(x),2), " (",round(sd(x),2), ")"))), 
        by=model
      ]
      # true values
      true.params = c("DGF", scen.coef[[scenario]]$value)
      names(true.params) = c("model", scen.coef[[scenario]]$variable)
      if (scenario == 1){
        dt.final = rbind(
          t(true.params),
          dt.table, 
          fill=T
        )
        dt.final$Scen = scenario
      } else {
        # JOINING 
        dt.tmp = rbind(
          t(true.params),
          dt.table, 
          fill=T
        )
        dt.tmp$Scen = scenario
        dt.final = rbind(dt.final, dt.tmp)
      } # B REPETITIONS
      print("HELLO I'm DONE :) ")
    } # CENSOR 
  } # N SAMPLE SIZE
  
  setcolorder(dt.final, c("Scen","model", mvars))
  fname = paste0("./documents/thesis_paper/sim_surv_table.tex")
  print(paste("write latex table to file",fname))
  print(
    xtable::xtable(
      x=dt.final, 
      caption=c("Survival Regression Parameter Estimates")
    ),
    file=fname,
    include.rownames=FALSE,
    floating=FALSE,
    # latex.environments = "flushleft",
    size="small",
  )
}


## RUN SIM ------
# sim.surv.cox.aft()

