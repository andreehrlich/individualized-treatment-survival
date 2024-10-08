##########################################
# ITR COMPARISON ON REAL DATA SETS 
##########################################
library(data.table) # data manipulation
library(radiant.data) # which.pmin()
library(survival) # methods + data sets 
library(kernlab)  # rbf kernel
library(pbapply)  # parallel
library(parallel)
# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
# source("./models/msowl/R/beta_version/msowl_new.R")
source("./models/msowl_fresh.R")


##############################
# MSOWL Provided EXAMPLE Data Set 
##############################
# spectrum_os = read.csv("./data/spectrum_os.csv")


data.example = read.csv("./models/msowl/data/example_data.csv")
my_formula = ms(t1 = t1, t2 = t2, s1 = s1, s2 = s2, a = A) ~ Z1 + Z2 
msowl_res = msowl(
  my_formula, 
  id = "id", 
  w = c(1,1,0), 
  tau = 3, 
  data.example, 
  lambda = 1,
  jackknife = FALSE, 
  trim = NULL
)
msowl_res
msowl_val = msowl.val(my_formula,"id", c(1,1,0),tau=3, data.example,msowl_res,fixed.rules=T,SE=T)
msowl_val
select.lambda(
  my_formula,
  id = "id",
  w = c(1, 1, 0), 
  tau = 3, 
  data = data.example,
  lambdas = c(0.01, 0.1, 0.5, 1, 5, 10, 50, 100)
)

#fit <- msowl(ms(t1 = t1, t2 = t2, s1 = s1, s2 = s2, a = A) ~ Z1 + Z2,
#             id = "id", w = c(0, 1, 0), tau = 3, data = dat,
#             lambda = 1)

#Value <- msowl.val(ms(t1 = t1, t2 = t2, s1 = s1, s2 = s2, a = A) ~ Z1 + Z2,
#                   id = "id", w = c(0, 1, 0), tau = 3, rule = fit, fixed.rules = TRUE, data = dat)



######################################################
## 2. retinopathy - preliminary analysis by Demiris
######################################################

## Time to vision loss
## Initially ignore within-patient correlation
##Dataset useful to illustrate random effects, i.e. two eyes per patient

fit <- survfit(Surv(futime, status) ~ trt, data = retinopathy)
plot(fit, col=1:2, xlab="Days from randomisation", ylab="Proportion with vision")

## adjust for eye

retinopathy$group[retinopathy$trt==0 & retinopathy$eye=="left"] <- 1
retinopathy$group[retinopathy$trt==1 & retinopathy$eye=="right"] <- 2
retinopathy$group[retinopathy$trt==0 & retinopathy$eye=="right"] <- 3
retinopathy$group[retinopathy$trt==1 & retinopathy$eye=="left"] <- 4
fit <- survfit(Surv(futime, status) ~ factor(group), data = retinopathy)
plot(fit, col=c(1,2,1,2), lty=c(1,1,2,2), xlab="Days from randomisation", ylab="Proportion surviving")
legend(5,0.3, c("right eye, not treated","right eye, treated", "left eye, not treated", "left eye, treated"), col=c(1,2,1,2), lty=c(1,1,2,2))
## Treatment appears to work better on right eye.

# ITR ANALYSIS
############################################

retinopathy

# hyper-parameters
w_single_state = c(1,0)
n = dim(retinopathy)[1]
num_cores = detectCores() - 1
ret.tau = 40

# convert data set for msowl analysis
ret_data = data.frame(
  id = 1:n,
  trt = as.factor(retinopathy$trt * 2 -1),
  status = as.factor(retinopathy$status),
  s1 = rep(1,n),
  s2 = 1 + retinopathy$status,
  t1 = rep(0,n),
  t2 = retinopathy$futime,
  age = retinopathy$age,
  laser = as.numeric(factor(retinopathy$laser))
  # = as.numeric(factor(retinopathy$type)),
  # eye = as.numeric(factor(retinopathy$eye))
)
ret_data
my_formula = ms(t1 = t1, t2 = t2, s1 = s1, s2 = s2, a = trt) ~ age + laser 


msowl_res = msowl(
    my_formula, 
    id = "id", 
    w = c(1,0), 
    tau = ret.tau, 
    ret_data, 
    lambda = 1,
    kernel="linear",
    # sigma=5,
    jackknife = FALSE, 
    trim = NULL,
    owl_method = "MSOWL",
    debug=T
)
  


# split data into test and train 
sample <- sample(
  c(TRUE, FALSE), 
  nrow(ret_data), 
  replace=TRUE, 
  prob=c(0.7,0.3)
)
train  <- ret_data[sample, ]
dim(train)
test   <- ret_data[!sample, ]
dim(test)

# tau = ??? 

library(ggplot2)

ggplot(
  ret_data, 
  aes(x = t2, fill = status, colour = status)) + 
  geom_histogram(alpha = 0.5, position = "identity") +
  ggtitle("Distribution of Survival Times by Event Status") 

# select.lambda(
#   my_formula,
#   id = "id", 
#   w = w_single_state, 
#   tau = ret.tau, 
#   data = train_data,
#   lambdas = c(0.01, 0.1, 0.5, 1, 5, 10, 50, 100)
# )
# compare itr estimators  
cresult = compare_itr_estimators(
  train, 
  test, 
  tau = ret.tau,
  my_formula = my_formula,
  num_cov = 4
)
cresult

hist(train$t2, breaks=100)

msowl_linear_fit <- msowl(
  formula=my_formula,
  id = "id",
  w = w_single_state,
  tau = 60,
  data = as.data.frame(train),
  kernel="linear",
  lambda = 1,
  jackknife = FALSE,
  trim = NULL
)



# ###############################
# 
# 
# 
# # formula for msowl
# my_formula = ms(t1 = t1, t2 = t2, s1 = s1, s2 = s2, a = A) ~ z1 + z2 + z3
# 
# ######################################
# # precalculate BIG-N Test set reward ?
# ######################################
# data.test.complete = FALSE
# if(n.test > 1000){
#   data.test.complete = precalculate_reward(
#     data=data.test,
#     n.test=n.test,
#     my_formula=my_formula,
#     tau=tau,
#     n.batch = n.batch
#   )
# }
# ######################################
# 
# ######################################
# # divide_and_conquer_value 
# ######################################
# # n_groups = n.test / 200
# # 
# # data.test.batches = split(
# #   as.data.frame(data.test),
# #   sample.int(
# #     n=n_groups,
# #     size=nrow(data.test),
# #     replace=TRUE
# #   )
# # )
# ####################################
# 
# ####################################
# if (debug) print("MSOWL")
# ####################################
# msowl_train <- msowl(
#   formula=my_formula,
#   id = "id",
#   w = c(1,0),
#   tau = tau,
#   data = as.data.frame(data.train),
#   # kernel="linear",
#   lambda = 1,
#   jackknife = FALSE,
#   trim = NULL,
#   debug=F
# )
# if (debug) print(msowl_train$Value)
# 
# # TODO MAKE PREDICTIONS FROM THE SVM !!! 
# # msowl_train$call
# # get treatment recommendations
# # form <- as.formula(paste("~", strsplit(as.character(msowl_train$call), "~", fixed = TRUE))[[3]])
# # Z <- model.matrix(form, dat)
# # f_n <- msowl_train$fit$Z%*%msowl_train$beta_opt
# # class(msowl_train)
# 
# #####
# if(is.list(data.test.complete)){
#   msowl_test = msowl.val(
#     formula = my_formula,
#     id = "id",
#     w = c(1,0),
#     tau = tau,
#     data = as.data.frame(data.test),
#     rule = msowl_train,
#     fixed.rules=T,
#     reward = data.test.complete
#   )
#   if (debug) print(msowl_test)
# } else {
#   #########
#   msowl_test = pbapply::pbsapply(
#     data.test.batches,
#     msowl.val,
#     formula=my_formula,
#     id="id",
#     w = c(1,0),
#     tau=tau,
#     rule=msowl_train,
#     fixed.rules=T
#   )
#   # print(mean(msowl_test))
# }
# 
# ####################################
# if (debug) print("ICO")
# ####################################
# ico_train = msowl(
#   formula=my_formula,
#   id = "id",
#   w = c(1,0),
#   tau = tau,
#   data = as.data.frame(data.train),
#   # kernel="linear",
#   lambda = 1,
#   jackknife = FALSE,
#   trim = NULL,
#   owl_method = "ICO",
#   debug=F
# )
# if(unlist(gregexpr('error', class(ico_train)))[1] !=-1){
#   # error
#   ico_train_val = NA
#   ico_test = NA
#   if (debug) print("ICO error!")
# } else{
#   ico_train_val = ico_train$Value
#   
#   #########
#   if(is.list(data.test.complete)){
#     
#     ico_test = msowl.val(
#       formula = my_formula,
#       id = "id",
#       w = c(1,0),
#       tau = tau,
#       rule = ico_train,
#       fixed.rules=FALSE,
#       data = as.data.frame(data.test),
#       reward = data.test.complete
#     )
#     if (debug) print(ico_test)
#   } else {
#     #########
#     ico_test = pbapply::pbsapply(
#       data.test.batches,
#       msowl.val,
#       formula=my_formula,
#       id="id",
#       w = c(1,0),
#       tau=tau,
#       rule=ico_train,
#       fixed.rules=F
#     )
#     #########
#   }
# }
# ####################################
# if (debug) print("COX")
# ####################################
# # specify cox regression formula 
# cox_form = formula=as.formula("Surv(time=t1, time2=t2, status.event) ~ z1 + z2 + z3 + A + A:z1 + A:z2 + A:z3")
# if (scenario %in% c(1,3)){
#   cox_form = as.formula("Surv(time=t1, time2=t2, status.event) ~ z1 + A + A:z1")
# } else if (scenario %in% c(2,4)){
#   cox_form = as.formula("Surv(time=t1, time2=t2, status.event) ~ A + A:z1 + A:z2")
# }
# 
# # event status = {1: Dead, 0: Censored}
# # states = {1: alive, 2: dead}
# data.train$status.event = as.integer(data.train$s2 == 2)
# 
# cox.train <- coxph(
#   formula=cox_form, 
#   data=data.train
# )
# # summary(cox.train.itr)
# 
# if(is.list(data.test.complete)){
#   # PRECALCULATED REWARD
#   
#   # EVAL TRAIN SET
#   cox.train.predictions <- coxph_itr_predict(
#     cox_fit = cox.train,
#     data = data.train,
#     tau=tau
#   )
#   cox.train.val <- msowl.val(
#     formula=my_formula,
#     id="id",
#     w = c(1,0),
#     tau=tau,
#     rule = cox.train.predictions,
#     fixed.rules=FALSE,
#     data=as.data.frame(data.train),
#     reward = data.test.complete
#   )
#   if (debug) print(cox.train.val)
#   
#   # EVAL TEST SET
#   cox.test.predictions <- coxph_itr_predict(
#     cox_fit = cox.train,
#     data = data.test,
#     tau=tau
#   )
#   cox.test.val <- msowl.val(
#     formula=my_formula,
#     id="id",
#     w = c(1,0),
#     tau=tau,
#     rule = cox.test.predictions,
#     fixed.rules=FALSE,
#     data=as.data.frame(data.test),
#     reward = data.test.complete
#   )
#   if (debug) print(cox.test.val)
#   ######################
# } else {
#   #  DIVIDE & CONQUER  
#   cox_dq_helper = function(
#     data,
#     cox_fit,
#     formula=my_formula,
#     id="id",
#     w = c(1,0),
#     tau=tau
#   ){
#     cox.test.predictions <- coxph_itr_predict(
#       cox_fit = cox_fit,
#       data = data,
#       tau=tau
#     )
#     
#     cval = msowl.val(
#       data=data,
#       formula=my_formula,
#       id="id",
#       w = c(1,0),
#       tau=tau,
#       rule=cox.test.predictions,
#       fixed.rules=F
#     )
#   }
#   
#   cox.test.val = pbapply::pbsapply(
#     X=data.test.batches,
#     FUN=cox_dq_helper,
#     cox_fit = cox.train,
#     formula=my_formula,
#     id="id",
#     w = c(1,0),
#     tau=tau
#   )
#   ############################################
# }
# ###################
# # RESULTS 
# ###################
# # train_result = c(
# #   msowl_linear_train=msowl_train$Value,
# #   ico_linear_train = ico_train_val
# # )
# # # train_result
# 
# test_result = c(
#   n=n.train,  
#   scen=scenario,
#   pct_cens =pct_censored_train,
#   b=b,
#   dur = Sys.time() - round_time)
# 
# if(is.list(data.test.complete)){
#   # precalculated reward
#   test_result = c(
#     test_result, 
#     PLUS =unname(msowl_test["V(1)"]),
#     MINUS=unname(msowl_test["V(-1)"]),
#     MSOWL=unname(msowl_test["V(dn)"]),
#     ICO = unname(ico_test),
#     COX=unname(cox.test.val["V(dn)"])
#   )
# } else{
#   # div & conquer
#   test_result = c(
#     test_result, 
#     PLUS.mu  = mean(msowl_test[2,]),
#     PLUS.sd  = sd(msowl_test[2,]),
#     MINUS.mu = mean(msowl_test[3,]),
#     MINUS.sd = sd(msowl_test[3,]),
#     MSOWL.mu = mean(msowl_test[1,]),
#     MSOWL.sd = sd(msowl_test[1,]),
#     ICO.mu   = mean(ico_test),
#     ICO.sd   = sd(ico_test),
#     COX.mu   = mean(cox.test.val),
#     COX.sd   = sd(cox.test.val)
#   )
# }
# if (debug) print(test_result)
# return(test_result)
# 
# 
