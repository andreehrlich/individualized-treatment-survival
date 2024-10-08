#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# ZhaoKosorok2015 ITR Simulation Experiment Replication
# MSc Statistics Thesis: Individualized Treatment Regimes (ITR)  2024
# Andre Ehrlich
# Supervisors: Demiris, Bakoyannis

# DESCRIPTION -----
# sim.kos():: Simulation Study:
# - Replication of Zhao Kosorok et al Doubly Robust 2015 
# - Compare performance of several ITR methods 
# - Vary treatment effect functions, censoring rates, sample sizes
# - Proportional hazard assumption vs not
#
# Survival Models 
# - Cox
# - AFT 
# - Exponential 
#
# ITR Models 
# - Standard Zero-Order "One-Size-Fits-All" Model
# - Backwards-Induction using Cox Regression
# - Outcome Weighted Learning (OWL)
# - - Multistate OWL (Bakoyannis2023)
# - - ICO OWL (ZhaoKosorok2015)
# - - Doubly-Robust OWL (ZhaoKosorok2015)

# LIBRARIES

# plot
library(ggplot2)
library(survival)
library(ggh4x)

# parallel 
library(pbapply)
pbo = pboptions(type="txt") # timer
library(parallel)

# library(future)
# library(future.apply)
# plan(multicore)

# models 
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# getwd()
source("./msowl_modified.R")
source("./sim_kos_plot.R")
source("./survival.regression.R")
library(radiant.data)
library(xtable)

# MAIN SIM LOOP ----------------------------------------------------------------
doOne = function(my_input) {
  round_time    = Sys.time()
  scenario      = my_input[[1]]
  censor        = my_input[[2]]
  n.train       = my_input[[3]]
  n.test        = my_input[[4]]
  n.batch       = my_input[[5]]
  tau           = my_input[[6]]
  b             = my_input[[7]]
  debug         = my_input[[8]]
  debug.tuning  = my_input[[9]]
  div.method    = my_input[[10]]
  kfolds        = my_input[[11]]
  lambda        = my_input[[12]]
  ico.lambda.in = my_input[[13]]
  SE.tune       = my_input[[14]]
  if (debug) {print("SIMULATION LOOP"); print(my_input)}
  
  formula = ms(t1 = t1, t2 = t2, s1 = s1, s2 = s2, a = A) ~ z1 + z2 + z3
  # TRAIN  -----------------------------------------------------------
  data.train = sim.kos(
    n           = n.train,
    scenario    = scenario,
    tau         = tau,
    censor.rate = censor
  )
  train.pct.censored = mean(!data.train$status)
  
  if (debug) print(data.train)
  if (debug) print(paste("train cens pct: ", train.pct.censored))
  # TEST  -----------------------------------------------------------
  data.test = sim.kos(
    n           = n.test,
    scenario    = scenario,
    tau         = tau,
    censor.rate = FALSE
  )
  test.pct.cens = mean(!data.test$status)
  if (debug) print(data.test)
  if (debug) print(paste("test cens pct: ",test.pct.cens))
  # PRECALC  -----------------------------------------------------------
  time.pre.calc = NA
  count.pre.calc = 0
  data.train.complete = NULL
  data.train.complete.ico = NULL
  data.test.complete = NULL
  # if (div.method == 1){
  batches = split(
    as.data.frame(data.test),
    sample.int(
      n=n.test/n.batch,
      size=nrow(data.test),
      replace=TRUE
    )
  )
  if (div.method == 2) {
    time.pre.calc = Sys.time()
    count.pre.calc = 0
    ## TRAIN --------
    if(n.train >= 3*n.batch){
      if (debug) print("PRECALC TRAINING DATA FOR EACH MODEL")
      count.pre.calc = count.pre.calc + 1
      data.train.complete = precalculate_reward(
        data = data.train,
        id = "id",
        w = c(1,0),
        formula = formula,
        tau = tau,
        n.batch = n.batch,
        owl_method = "MSOWL"
      )
      count.pre.calc = count.pre.calc + 1
      data.train.complete.ico = precalculate_reward(
        id = "id",
        w = c(1,0),
        data = data.train,
        formula = formula,
        tau = tau,
        n.batch = n.batch,
        owl_method = "ICO"
      )
    }
    ## TEST --------
    if(n.test >= 3*n.batch){
      if (debug) print("PRECALC TEST DATA")
      count.pre.calc = count.pre.calc + 1
      data.test.complete = precalculate_reward(
        data=data.test,
        id = "id",
        w = c(1,0),
        # n=n.test,
        formula=formula,
        tau=tau,
        n.batch = n.batch
      )
    }
    time.pre.calc = Sys.time() - time.pre.calc
    if (debug) print(paste("PRE CALC: ", time.pre.calc, count.pre.calc))
  }
  ## MODEL: MSOWL --------------------------------------------------
  if (debug) print("MSOWL")
  ### LAMBDA TUNING ----
  if (lambda == FALSE) {
    if (debug) print("MSOWL TUNING KFOLD")
    lambda.tuning = select.lambda.kfold(
      formula = formula,
      id = "id",
      w = c(1,0),
      tau = tau,
      data = as.data.frame(data.train),
      lambdas = lambda.vec,
      owl_method = "MSOWL",
      k = kfolds,
      round_vals = NULL,
      debug = debug.tuning,
      SE = SE.tune
    )
    if (debug) print(lambda.tuning)
    msowl.lambda = dim(data.train)[[1]]^(-1/2) * lambda.tuning$lambda
  } else if (lambda == -1) {
    if (debug) print("MSOWL TUNING LOOCV ")
    lambda.tuning = select.lambda(
      formula, 
      id = "id", 
      w = c(1,0), 
      tau=tau, 
      as.data.frame(data.train), 
      lambdas = lambda.vec,
      owl_method = "MSOWL",
      debug = debug.tuning
    )
    if (debug) print(lambda.tuning)
    msowl.lambda = dim(data.train)[[1]]^(-1/2) * lambda.tuning$lambda
  } else {
    msowl.lambda  = dim(data.train)[[1]]^(-1/2) * lambda
    if (debug) print(paste0("lambda manually selected: ", lambda))
  }
  ### TRAIN ---------
  msowl_train <- msowl(
    formula = formula,
    id = "id",
    w = c(1,0),
    tau = tau,
    data = as.data.frame(data.train),
    lambda = msowl.lambda, 
    jackknife = FALSE,
    trim = NULL,
    reward = data.train.complete,
    debug=F
  )
  if (debug) {print("msowl_train"); print(str(msowl_train))}
  if(unlist(gregexpr('error', class(msowl_train)))[1] !=-1){ 
    # ERROR
    itr.msowl.val = list('V(dn)'=NA, "V(1)"=NA, "V(-1)"=NA)
    time.msowl.val = NA
  } else{
    ### TEST ------------
    if (debug) print(msowl_train$Value)
    time.msowl.val = Sys.time()
    if (div.method == 1) {
      res = pbapply::pbsapply(
        batches,
        msowl.val,
        formula=mformula,
        id="id",
        w = c(1,0),
        tau=tau,
        rule=msowl_train,
        fixed.rules=T
      )
      itr.msowl.val = c(
        "V(dn)" = mean(res[1,]),
        "V(-1)" = mean(res[2,]),
        "V(1)" = mean(res[3,])
      )
    } else if (div.method == 2) {
      itr.msowl.val = msowl.val(
        data = as.data.frame(data.test),
        rule = msowl_train,
        formula = formula,
        id = "id",
        w = c(1,0),
        tau = tau,
        fixed.rules=T,
        reward=data.test.complete
      )
    }
    time.msowl.val = Sys.time() - time.msowl.val
    if (debug) {
      # print(itr.msowl.val)
      # print(time.msowl.val)
      print(paste(c("MSOWL", itr.msowl.val, time.msowl.val)))
    }
  }
  ## MODEL: ICO --------------------------------------------------
  ### LAMBDA TUNING -----------
  if (ico.lambda.in == FALSE) {
    if (debug) print("tuning ICO")
    ico.lambda.tuning = select.lambda.kfold(
      formula = formula,
      id = "id",
      w = c(1,0),
      tau = tau,
      data = as.data.frame(data.train),
      lambdas = lambda.vec,
      owl_method = "ICO",
      k = kfolds,
      round_vals = NULL,
      debug = debug.tuning,
      SE = SE.tune
    )
    if (debug) print(ico.lambda.tuning)
    ico.lambda  = dim(data.train)[[1]]^(-1/2) * ico.lambda.tuning$lambda
  } else if (ico.lambda.in == -1) {
    ico.lambda.tuning = select.lambda(
      formula, 
      id = "id", 
      w = c(1,0), 
      tau=tau, 
      as.data.frame(data.train), 
      lambdas = lambda.vec,
      owl_method = "ICO",
      debug = debug.tuning
    )
    ico.lambda  = dim(data.train)[[1]]^(-1/2) * ico.lambda.tuning$lambda
  } else {
    ico.lambda  = dim(data.train)[[1]]^(-1/2) * ico.lambda.in
  }
  
  ###  TRAIN -------
  if (debug) print("ICO")
  ico_train = msowl(
    formula=formula,
    id = "id",
    w = c(1,0),
    tau = tau,
    data = as.data.frame(data.train),
    lambda = ico.lambda,
    jackknife = FALSE,
    trim = NULL,
    owl_method = "ICO",
    reward = data.train.complete.ico,
    debug=F
  )
  if(unlist(gregexpr('error', class(ico_train)))[1] !=-1){
    # ERROR
    ico_train_val = NA
    itr.ico.val = NA
    time.ico.val = 0
    if (debug) print("ICO error!")
  } else {
    ### TEST -------
    time.ico.val = Sys.time() 
    ico_train_val = ico_train$Value
    if (div.method == 1) {
      res = pbapply::pbsapply(
        batches,
        msowl.val,
        formula=formula,
        id="id",
        w = c(1,0),
        tau=tau,
        rule=ico_train,
        fixed.rules=F
      )
      itr.ico.val = c("V(dn)" = mean(res))
    } else if (div.method == 2) {
      itr.ico.val = msowl.val(
        formula = formula,
        id = "id",
        w = c(1,0),
        tau = tau,
        rule = ico_train,
        fixed.rules=FALSE,
        data = as.data.frame(data.test),
        reward = data.test.complete
      )
    }
    time.ico.val = Sys.time() - time.ico.val
    if (debug) print(paste("ICO", itr.ico.val, time.ico.val))
  }
  ## MODEL: COX -----------------------------------------------------------
  ### TRAIN --------------
  if (debug) print("COX")
  cox.formula = as.formula("Surv(time=t1, time2=t2, status) ~ z1 + z2 + z3 + A + A:z1 + A:z2 + A:z3")
  # data.train$status = as.integer(data.train$s2 == 2)
  cox.train <- coxph(
    formula=cox.formula, 
    data=data.train
  )
  
  summary(cox.train)
  if (debug) summary(cox.train)
  ###  TEST  --------------
  time.cox.val = Sys.time()
  if (div.method == 1) {
    itr.cox.val = pbapply::pbsapply(
      X=batches,
      FUN=itr.cox.divide.conquer.value,
      cox_fit=cox.train,
      formula=formula,
      id="id",
      w = c(1,0),
      tau=tau
    )
    itr.cox.val = c("V(dn)" = mean(itr.cox.val))
  } else if (div.method == 2) {
    cox.test.predictions <- coxph_itr_predict(
      cox_fit = cox.train,
      data = data.test,
      tau=tau
    )
    itr.cox.val <- msowl.val(
      formula=formula,
      id="id",
      w = c(1,0),
      tau=tau,
      rule = cox.test.predictions,
      fixed.rules=FALSE,
      data=as.data.frame(data.test),
      reward = data.test.complete
    )
  }
  time.cox.val = Sys.time() - time.cox.val
  if (debug) print(paste("COX", itr.cox.val, time.cox.val))
  ## AGG RESULTS --------------------------------------------------
  test_result = c(
    # simulation params
    n           = n.train,  
    t.test      = n.test,
    scen        = scenario,
    cens_param  = censor,
    pct.cens    = round(train.pct.censored,2),
    b           = b,
    dur         = round(Sys.time() - round_time,2),
    COMP_METHOD = div.method,
    lambda      = msowl.lambda,
    ico.lambda  = ico.lambda,
    # Values
    PLUS        = unname(itr.msowl.val["V(1)"]),
    MINUS       = unname(itr.msowl.val["V(-1)"]),
    MSOWL       = unname(itr.msowl.val["V(dn)"]),
    ICO         = unname(itr.ico.val["V(dn)"]),
    COX         = unname(itr.cox.val["V(dn)"]),
    # value comp durations 
    MSOWL.dur   = time.msowl.val,
    ICO.dur     = time.ico.val,
    COX.dur     = time.cox.val
    # PRE.dur   = time.pre.calc,
    # PRE.cnt   = count.pre.calc,
    # MSOWL.LAM.max = lambda.max,
    # MSOWL.LAM.dur = lambda.tuning.timer,
    # ICO.LAM.max = ico.lambda.max,
    # ICO.LAM.dur = ico.lambda.tuning.timer
  )
  if (debug) print(test_result)
  
  # if (test_result['MSOWL'] < Vest_result['MINUS'] | test_result['MSOWL'] < test_result['PLUS']){
  #   print("ERROR: MSOWL UNDERPERFORMS ZOM!!!")
  # }
  return(test_result)
}


# MC SIM ----------------------------------------------------------------------
run.sim.kos <- function(
  B = 100,
  n.test = 10000,
  n.batch = 200,
  n.vec = c(100, 200, 400),
  scenario.vec = c(1, 2, 3, 4),
  censor.vec = c(-1, -0.5, -0.2),
  keyword = "KOS_SIM",
  kfolds = 5,
  div.method   = 2,
  lambda = 0.01,
  ico.lambda = 0.01,
  SE.tune = F,
  test_one = T
) {
  
  
  cat( "KOS SIMULATION", "B", B, "n.test", n.test, "n.batch", n.batch, 
       "scenario.vec", scenario.vec, "censor.vec", censor.vec, 
       "n.vec", n.vec, "keyword", keyword, "kfolds", kfolds, 
       "div.method", div.method, "lambda", lambda, "ico.lambda", 
       ico.lambda, "SE.tune", SE.tune, "test_one", test_one)
  
  foldername = paste0("./sim_results/", keyword)
  dir.create(foldername, showWarnings = FALSE)
  print(paste0("Results to be saved in ", foldername))
  
  ## PARALLEL  -----
  switch(
    Sys.info()[['sysname']],
    Windows= {
      # remote compute
      print("I'm a Windows PC. Use all cores")
      num_cores = detectCores()
      print(paste("num_cores", num_cores))
    },
    Darwin = {
      # personal computer
      print("I'm a Mac. Use all cores but one.")
      # source("./models/msowl/R/beta_version/msowl_new.R")
      num_cores = detectCores()-1
      print(paste("num_cores", num_cores))
    }
  )
  # CONFIGURATION (FOR BOTH WINDOWS & MAC/LINUX)
  cl = makeCluster(
    num_cores,
    outfile="" # paste0("./sim_results/", keyword,"/sim_log.txt")
  )
  clusterExport(
    cl,
    varlist = c(
      # survival analysis
      "coxph",
      "Surv",
      "survfit",
      # my functions 
      "sim.kos",
      "sim.surv.regr",
      "H.inv.kos",
      "select.lambda.kfold",
      "select.lambda",
      "coxph_itr_predict",
      # msowl
      "precalculate_reward",
      "calculate_reward",
      "ms",
      "zeta_t",
      "basehaz",
      "bind_rows",
      "msowl",
      "msowl.val",
      "wsvm_solve",
      "ipop",
      "primal",
      "lambda.vec",
      # helpers
      "copy",
      "which.pmax",
      "which.pmin",
      "pbsapply"
    )
  )

  for (n.train in n.vec){
    for (censor in censor.vec) {
      for(scenario in scenario.vec){
        # FNAME ------
        c_fname = paste0(foldername,"/c",censor,"_s",scenario,"_n", n.train,".csv")
        print(c_fname)
        if (file.exists(c_fname)){
          print("this file already exists. skipping.")
          next
        }
        ### INPUT ARGS ---------
        my_inputs  = data.frame(
          scenario     = rep(scenario,B),
          censor       = rep(censor,B),
          n.train      = rep(n.train,B),
          n.test       = rep(n.test,B),
          n.batch      = rep(n.batch,B),
          tau          = rep(tau.vec[scenario],B),
          b            = 1:B,
          debug        = rep(F,B),
          debug.tuning = rep(F,B),
          div.method   = rep(div.method,B),
          kfolds       = rep(kfolds,B),
          lambda       = rep(lambda,B),
          ico.lambda   = rep(ico.lambda,B),
          SE.tune      = rep(SE.tune,B)
        )
        
        if (test_one) {
          print("TEST ONCE")
          # print(my_inputs[1,])
          ind = sample.int(nrow(my_inputs), 1)
          print(paste("ind",ind))
          test_input = copy(my_inputs[ind,])
          test_input$debug = TRUE
          # test_input$scenario = 4
          # test_input$tau = tau.vec[4]
          print(test_input)
          my_res = doOne(test_input)
          print(my_res)
        }
        
        ### PBAPPLY -----
        print("tail(my_inputs)")
        print(tail(my_inputs))
        print(dim(my_inputs))
        print(paste("parrallel start time:", Sys.time()))
        cres = pbapply(
          X = my_inputs, 
          FUN=doOne, 
          MARGIN=1,
          cl= cl
        )
        ### WRITE TO FILE -------
        print("write to file.")
        write.table(
          t(cres),
          c_fname,
          append = FALSE,
          row.names=FALSE,
          col.names=TRUE,
          sep=","
        )
      } # scen
    } # cens
  } # n.train
  stopCluster(cl)
}

# RUN SCRIPT ---------
## RUN FULL SIM ---------
if (length(args) > 0) {
  if (args[1] == "sim") {
    keyword="run16"
    n.vec = c(100)    # ,200 ,400)
    censor.vec = c(F) # ,-1,  -0.5)
    scenario.vec = c(1,2,3,4)
    run.sim.kos(
      keyword      =  keyword,
      B            =  100,
      n.vec        =  n.vec,
      n.test       =  10000,
      censor.vec   =  censor.vec,
      scenario.vec =  scenario.vec,
      kfolds       =  5,
      div.method   =  2,
      lambda       = -1,
      ico.lambda   = -1,
      SE.tune      =  F,
      test_one     =  T
    )
    ## PLOT RESULTS ------
    dt.kos = setDT(
      load.files.from.dir(
        keyword      = keyword,
        n.train.vec  = n.vec,
        cens.vec     = censor.vec,
        scen.vec     = scenario.vec
      )
    )
    plot.sim.big.boxplots(dt.kos, keyword)
  } else if (args[1] == "test") {
    # Test Once -----------
    my_input = c(
      scenario = 1,
      censor = F,
      n.train = 100,
      n.test = 10000,
      n.batch = 200,
      tau = tau.vec[1],
      b = 1,
      debug = T,
      debug.tuning = T,
      div.method = 2,
      kfolds = 5,
      lambda = F,
      ico.lambda.in = F,
      SE.tune   = T
    )
    res_debug = doOne(
      my_input = my_input
    )
  }
}