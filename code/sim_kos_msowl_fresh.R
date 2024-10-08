# Kosorok ITR Simulation Replication
# MSc Statistics Thesis: Individualized Treatment Regimes (ITR)  2024
###############################################################################
# Author: Andre Ehrlich
# Supervisors: Bakoyannis, Demiris 
###############################################################################
# Simulation Study:
# - Compare performance of several ITR methods 
# - Vary treatment effect functions, censoring rates, sample sizes
# - Proportional hazard assumption vs not
# ITR Models under comparison
# - Standard Zero-Order "One-Size-Fits-All" Model
# - Backwards-Induction using Cox Regression
# - Outcome Weighted Learning (OWL)
# - - Multistate OWL (Bakoyannis2023)
# - - ICO OWL (ZhaoKosorok2015)
# - - Doubly-Robust OWL (ZhaoKosorok2015)
###############################################################################
library(pbapply)      # parallel
library(parallel)
library(radiant.data) # which.pmin()

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
# source("./msowl_modified.R")
source("../msowl_fresh/R/beta_version/msowl_new.R")

# preferences between windows / mac 
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

###########################################
# simulation hyperparameters
###########################################
# num_cores = detectCores()-1
B = 1000
scenario.vec = c(1,2,3,4)
tau.vec = c(1.5, 2, 2, 2.5)
n.vec = c(400,100,200,800,1600)
n.test = 10000
n.batch = 200
censor.vec = c(-1, -0.5, -0.2)
#######################################################################
# Simulation: Replication of Zhao Kosorok et al Doubly Robust 2015 
#######################################################################
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
  time.vec.censoring = rexp(n,rate= exp(censor.rate))
  if (censor.rate == FALSE) time.vec.censoring = rep(10000,n) 
  # SURVIVAL
  if (scenario == 1){
    time.vec.survival = sim.cox(
      n = n, 
      H_inv = H_inv_kos, 
      psi = 0.6*z1 - 0.8*z2 + (0.6 - 0.4*z1 - 0.2*z2 - 0.4*z3)*A
    )
  } else if (scenario == 2){
    time.vec.survival = sim.cox(
      n=n, 
      H_inv=H_inv_kos, 
      psi = -1.5*z1 + 0.5*z2 + (1 - 0.7*z1 - 1.2*z2)*A
    )
  } else if (scenario == 3){
    error.term = rnorm(n, 0,1)
    linear.predictor = -0.5 - 0.8*z1 + 0.7*z2 + 0.2*z3 +(0.6-0.4*z1 -0.1*z2-0.4*z3)*A
    time.vec.survival = exp(linear.predictor + error.term)
    
  } else if (scenario == 4){
    error.term = rnorm(n, 0,1)
    linear.predictor = -0.2 - 0.5*z1 + 0.5*z2 + 0.3*z3 +(0.5-0.1*z1 -0.6*z2+0.1*z3)*A
    time.vec.survival = exp(linear.predictor + error.term)
  }
  # OBSERVED INFORMATION 
  time.vec.observed = pmin(time.vec.survival, time.vec.censoring)
  status.event = as.integer(time.vec.survival <= time.vec.censoring)
  state.terminal = status.event + 1   # state 1 = alive; state 2 = dead
  records = data.frame(
    id = 1:n,
    z1 = z1,
    z2 = z2,
    z3 = z3,
    A  = A,
    s1 = 1,
    s2 = state.terminal,
    t1 = 0,
    t2 = time.vec.observed
    # d = status.event
  )
  records
  pct_censored = mean(!status.event)
  return(list(records, pct_censored))
}

###########################
# simulation loop   
###########################
doOne = function(
    my_input
) {
  round_time = Sys.time()
  # PARAMATERS 
  scenario = my_input[[1]]
  censor   = my_input[[2]]
  n.train  = my_input[[3]]
  n.test   = my_input[[4]]
  n.batch  = my_input[[5]]
  tau      = my_input[[6]]
  b        = my_input[[7]]
  debug    = my_input[[8]]
  if (debug) {print("DEBUG MODE"); print(c(my_input))}
  
  my_formula = ms(t1 = t1, t2 = t2, s1 = s1, s2 = s2, a = A) ~ z1 + z2 + z3
  svm.lambda = n.train^(-1/2) # regularity conditions
  # TRAIN
  res = sim.kos(
    n = n.train,
    scenario = scenario,
    tau = tau,
    censor.rate = censor
  )
  data.train = res[[1]]
  pct_censored_train = res[[2]]
  if (debug) print(paste("train cens pct: ",pct_censored_train))
  
  # TEST
  res = sim.kos(
    n = n.test,
    scenario = scenario,
    tau = tau,
    censor.rate = FALSE
  )
  data.test = res[[1]]
  pct_censored = res[[2]]
  if (debug) print(paste("test cens pct: ",pct_censored))
  ############################################
  # PRECOMPUTE REWARD FOR "LARGE" DATA SETS
  ############################################
  data.train.complete = NULL
  data.train.complete.ico = NULL
  data.test.complete = FALSE
  # if(n.train > 3*n.batch){
  #   data.train.complete = precalculate_reward(
  #     data=data.train,
  #     n=n.train,
  #     my_formula=my_formula,
  #     tau=tau,
  #     n.batch = n.batch
  #   )
  #   # print("data.train.complete")
  #   # print(data.train.complete)
  # }
  # if(n.train > 3*n.batch){
  #   data.train.complete.ico = precalculate_reward(
  #     data=data.train,
  #     n=n.train,
  #     my_formula=my_formula,
  #     tau=tau,
  #     n.batch = n.batch,
  #     owl_method="ICO"
  #   ) 
  # } 
  # # TEST 
  # if(n.test > 3*n.batch){
  #   data.test.complete = precalculate_reward(
  #     data=data.test,
  #     n=n.test,
  #     my_formula=my_formula,
  #     tau=tau,
  #     n.batch = n.batch
  #   )
  # }
  ####################################
  if (debug) print("MSOWL")
  msowl_train <- msowl(
    formula=my_formula,
    id = "id",
    w = c(1,0),
    tau = tau,
    data = as.data.frame(data.train),
    # kernel="linear",
    lambda = svm.lambda,
    jackknife = FALSE,
    trim = NULL
    # reward = data.train.complete,
    # debug=F
  )
  if (debug) print(paste0("msowl_train: ",msowl_train$Value))
  if(unlist(gregexpr('error', class(msowl_train)))[1] !=-1){ # ERROR
    msowl_test = list('V(dn)'=NA, "V(1)"=NA, "V(-1)"=NA)
  } else{
    msowl_test = msowl.val(
      formula = my_formula,
      id = "id",
      w = c(1,0),
      tau = tau,
      data = as.data.frame(data.test),
      rule = msowl_train,
      fixed.rules=T
      # reward = data.test.complete
    )
    if (debug) print(msowl_test)
  }
  ####################################
  if (debug) print("ICO")
  ico_train = msowl(
    formula=my_formula,
    id = "id",
    w = c(1,0),
    tau = tau,
    data = as.data.frame(data.train),
    # kernel="linear",
    lambda = svm.lambda,
    jackknife = FALSE,
    trim = NULL,
    owl_method = "ICO"
    # reward = data.train.complete.ico,
    # debug=F
  )
  if(unlist(gregexpr('error', class(ico_train)))[1] !=-1){
    # ERROR
    ico_test = NA
    if (debug) print("ICO error!")
  } else{
    ico_train_val = ico_train$Value
    if (debug) print(paste0("ico_train: ",ico_train_val))
    ico_test = msowl.val(
      formula = my_formula,
      id = "id",
      w = c(1,0),
      tau = tau,
      rule = ico_train,
      fixed.rules=FALSE,
      data = as.data.frame(data.test)
      # reward = data.test.complete
    )
    if (debug) print(paste0("ice_test: ",ico_test))
  }
  ####################################
  if (debug) print("COX")
  cox.formula = as.formula("Surv(time=t1, time2=t2, status.event) ~ z1 + z2 + z3 + A + A:z1 + A:z2 + A:z3")
  data.train$status.event = as.integer(data.train$s2 == 2)
  cox.train <- coxph(
    formula=cox.formula, 
    data=data.train
  )
  if (debug) summary(cox.train)
  # EVAL TEST SET
  cox.test.predictions <- coxph_itr_predict(
    cox_fit = cox.train,
    data = data.test,
    tau=tau
  )
  cox.test.val <- msowl.val(
    formula=my_formula,
    id="id",
    w = c(1,0),
    tau=tau,
    rule = cox.test.predictions,
    fixed.rules=FALSE,
    data=as.data.frame(data.test)
    # reward = data.test.complete
  )
  if (debug) print(cox.test.val)
  ####################################
  if (debug) print("COX-INTX")
  cox.formula = as.formula("Surv(time=t1, time2=t2, status.event) ~  A:z1 + A:z2 + A:z3")
  # data.train$status.event = as.integer(data.train$s2 == 2)
  cox.intx <- coxph(
    formula=cox.formula, 
    data=data.train
  )
  if (debug) summary(cox.intx)
  # EVAL TEST SET
  cox.intx.test.predictions <- coxph_itr_predict(
    cox_fit = cox.intx,
    data = data.test,
    tau=tau
  )
  cox.intx.test.val <- msowl.val(
    formula=my_formula,
    id="id",
    w = c(1,0),
    tau=tau,
    rule = cox.intx.test.predictions,
    fixed.rules=FALSE,
    data=as.data.frame(data.test)
    # reward = data.test.complete
  )
  if (debug) print(cox.intx.test.val)
  
  
  ############################################
  test_result = c(
    n     = n.train,  
    scen  = scenario,
    cens_param = censor,
    pct_cens =round(pct_censored_train,2),
    b     = b,
    dur   = round(Sys.time() - round_time,2),
    PLUS  = unname(msowl_test["V(1)"]),
    MINUS = unname(msowl_test["V(-1)"]),
    MSOWL = unname(msowl_test["V(dn)"]),
    ICO   = unname(ico_test),
    COX   = unname(cox.test.val["V(dn)"]),
    COX.INTX = unname(cox.intx.test.val["V(dn)"]) 
  )
  if (debug) print(test_result)
  return(test_result)
}



###########################################
# DEBUGGING
my_input = c(
  scenario=1,
  censor = -1.4,
  n.train=400,
  n.test=10000,
  n.batch=200,
  tau=tau.vec[1],
  b=1,
  debug=T
)
# run func once 
res_debug = doOne(
  my_input = my_input
)
print(res_debug)

# see variation in lambda in variation from randomly generated data set 
# # text select.lambda
# din = c()
# for (b in 1:100){
#   din = rbind(din, my_input)
# }
# cres = pbapply(
#   X=din,
#   FUN=doOne,
#   MARGIN=1,
#   cl=detectCores()
# )
# 
# ccres = c()
# for (b in 1:100){
#   ccres = c(ccres, cres[1]$my_input$lambda)
# }


# WINDOWS MACHINE CONFIGURATION 
# cl - A cluster object created by makeCluster, or an integer to indicate number of child-processes 
# *** integer values are ignored on Windows) for parallel evaluations ***
cl = makeCluster(
  num_cores,
  outfile=""
)

clusterExport(
  cl,
  varlist = c(
    "sim.kos",
    "sim.cox",
    "H_inv_kos",
    "precalculate_reward",
    "calculate_reward",
    "coxph",
    "Surv",
    "survfit",
    "ms",
    "zeta_t",
    "basehaz",
    "bind_rows",
    "msowl",
    "msowl.val",
    "wsvm_solve",
    "ipop",
    "primal",
    "coxph_itr_predict",
    "which.pmax",
    "select.lambda",
    "select.lambda.new"
  )
)

# save to folder.
foldername = "./sim_results/alpha"
dir.create(foldername, showWarnings = FALSE)
saved_files = c()
for (n.train in n.vec){
  for (censor in censor.vec) {
    for(scenario in rev(scenario.vec)){
      
      # parellel in sequential chunks to save results intermittently :) 
      c_fname = paste0(foldername,"/c",censor,"_s",scenario,"_n", n.train, ".csv")
      print(c_fname)
      
      if (file.exists(c_fname)){
        print("this file already exists. skipping.")
        next
      }
      saved_files = c(saved_files,c_fname)
      
      ###############################
      
      my_inputs = data.frame(
        scenario=rep(scenario,B),
        censor=rep(censor,B),
        n.train=rep(n.train,B),
        n.test=rep(n.test,B),
        n.batch=rep(n.batch,B),
        tau=rep(tau.vec[scenario],B),
        b=1:B,
        debug=rep(F,B)
      )
      cres = pbapply(
        X = my_inputs, 
        FUN=doOne, 
        MARGIN=1,
        cl= cl
      )
      
      ############################
      print("write to file.")
      write.table(
        t(cres),
        c_fname,
        append = FALSE,
        row.names=FALSE,
        col.names=TRUE,
        sep=","
      )
      ############################
    }
  }
}
stopCluster(cl)

# combine saved_files into single csv 

# using precalculated reward
