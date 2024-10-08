# sim_bak.R script does the following:
# 1) ADAPT SIMULATIONS FROM BAKOYANNIS 2023
# 2) Demonstrate DIVIDE & CONQUER for MSOWL Value Estimation.

# TODO Replace with original msowl code
source("./msowl_modified.R")
library(pbapply)      # parallel
library(parallel)
library(data.table) 
library(radiant.data) # which.pmin()
library(survival) 
library(kernlab)  

w_single_state = c(1,0)
  
simulate_single_state_data = function(
    n = 500,                                # study population
    tau = 3,                                # total observation length
    censoring_rate = -1.6,                  # right-censoring rate
    chosen_scenario = 1,                    # true decision rule scenario  
    prob_treat = 0.5                        # propensity
){
  A = rbinom(n=n, 1, prob = prob_treat) * 2 - 1
  z1 = runif(n,-1,1)
  z2 = runif(n,-1,1)
  # survival transition intensities adapted from Bakoyannis2023 MSOWL paper
  f_star = function(z1, z2, scenario) {
    switch(as.character(scenario),
     `1` = z1 + z2,
     `2` = 2 * (z1 - z2),
     `3` = 1 + z2 - exp(-z1),
     `4` = 2 * log(2 - z1 - z2) - 1.4
  )}
  transition_intensity = exp(-0.5 * z1 + 0.5 * z2 + A * f_star(z1, z2, chosen_scenario))
  survival_times = rexp(n,transition_intensity) 
  censoring_time = rexp(n,rate= exp(censoring_rate))
  survival_observed = pmin(survival_times,censoring_time)
  end_state = which.pmin(censoring_time, survival_times)
  records = data.frame(
    id = 1:n,
    z1 = z1,
    z2 = z2,
    A  = A,
    s1 = rep(1,n),
    t1 = rep(0,n),
    s2 = end_state,
    t2 = survival_observed
  )
  records
  return(records)
}

demonstrate.divide.conquer = function(
  B = 10,
  scenario = 1,
  censoring_rate = -0.4,
  n = 400,
  row_count = 0,
  prob_treat = 0.5,
  n_mc_integration = 10000,
  num_cores=detectCores() 
) {
  # regression formula specifying which variables are covariates, treatment, states, and time-points
  ms_formula = ms(t1 = t1, t2 = t2, s1 = s1, s2 = s2, a = A) ~ z1 + z2
  # TRAIN DATA 
  train_data = simulate_single_state_data(
    n = n,                            
    tau = tau,                        
    censoring_rate = censoring_rate,  
    chosen_scenario = scenario,       
    prob_treat = 0.5
  )
  # TEST DATA  
  test_data = simulate_single_state_data(
    n = n_mc_integration,
    tau = tau,
    censoring_rate = censoring_rate,
    chosen_scenario = scenario,
    prob_treat = 0.5
  )
  str(test_data)
  ##############################################
  print("ITR Method = MSOWL with Linear Kernel")
  # learn rule on train set
  msowl_linear_fit <- msowl(
    formula=ms_formula,
    id = "id",
    w = w_single_state,
    tau = tau,
    data = as.data.frame(train_data),
    kernel="linear",
    lambda = 1,
    jackknife = FALSE,
    trim = NULL
  )
  print("first we will do div&conq approach, since it only takes 45 seconds.")
  print("APPROACH #1) DIVIDE & CONQUER")
  ##########################################################################
  # Split the test data into chunks of 500, and process chunks sequentially 
  ##########################################################################
  group_size = 500
  n_groups = n_mc_integration / group_size 
  test_subsets = split(as.data.frame(test_data), sample(1:n_groups,nrow(test_data), replace=TRUE))
  # n=10k, batchsize=500, total calculation = ~45 seconds
  cres = pbapply::pbsapply(
    test_subsets,
    msowl.val,
    formula = ms_formula,
    id = "id",
    w = w_single_state,
    tau = tau,
    rule = msowl_linear_fit,
    fixed.rules=FALSE,
    cl= 1
  )
  mean(cres)
  #######################################
  print("Then run value func on complete 10k data set. It will take hours.")
  print("APPROACH #2) Naive Estimate Value of rule on big  test set")
  #######################################
  # inside of msowl.val calculate_reward, can get a time estimate 
  # of the zeta_t func application by using pbapply:pbsapply 
  # instead of sapply(tt, zeta_t)
  # n=10k ==> estimated ~2hrs
  
  timer = Sys.time()
  print(timer)
  msowl_linear_val = msowl.val(
    formula = ms_formula,
    id = "id",
    w = w_single_state,
    tau = tau,
    rule = msowl_linear_fit,
    fixed.rules=TRUE,
    data = as.data.frame(test_data)
  )
  msowl_linear_val
  timer = Sys.time() - timer
  print(timer)
}
# RUN SIMULATION FUNCTION
demonstrate.divide.conquer()
