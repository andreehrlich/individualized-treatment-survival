## Bakoyiannis 2023 Optimal ITR Example
# https://github.com/gbakoyannis/msowl

source("~/aueb/mowl/R/itr.R")
source("~/aueb/mowl/R/V_d.R")
source("~/aueb/mowl/R/wsvm_solve.R")
source("~/aueb/mowl/R/zeta_h_t.R")
source("~/aueb/mowl/R/zeta_t.R")


options(error=traceback) 

library(survival)




##########################################################################
# NOTE: This code exploring the data set was provided by Professor Demiris 
##########################################################################
## 2. Retinopathy

## Time to vision loss
## Initially ignore within-patient correlation
##Dataset useful to illustrate random effects, i.e. two eyes per patient

head(retinopathy) # matched eyes comparison
help(retinopathy)


########################################################################
# Brief Description from Demiris 
########################################################################
# A trial of laser coagulation as a treatment to delay diabetic retinopathy.
# 394 eyes, 197 patients were a 50% random sample of the patients with 
# "high-risk" diabetic retinopathy. One eye randomized to laser treatment and 
# the other eye received no treatment. Time from initiation of treatment to 
# the time when visual acuity dropped below 5/200 two visits in a row. 
# Built-in lag time of ~6 months (visits were every 3 months). 
# Survival times are time to vision loss in months,less 6.5 months lag. 
# Censoring was caused by death, dropout, or end of the study.

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

retinopathy


ret_data = NULL #  retinopathy

ret_data = data.frame(
  id = retinopathy$id, 
  A = retinopathy$trt,
  t1 = 0,
  t2 = retinopathy$futime,
)
retinopathy

  
  
fit <- survfit(Surv(futime, status) ~ trt, data = retinopathy)
plot(fit, col=1:2, xlab="Days from randomisation", ylab="Proportion with vision")



############################################

# The artificial dataset `example_data.csv` (included in this repository) contains observations from an illness-death process without recovery. The dataset can be obtained as follows

library(foreign)
# data <- read.csv("~/mowl/data/example_data.csv")
head(data)
#   id t1        t2         s1 s2  Z1          Z2          A
# 1  1 0.0000000 0.01380424  1  3 -0.80030228 -0.02345331  1
# 2  2 0.0000000 0.25340078  1  3 -0.77119721 -0.34741893  1
# 3  3 0.0000000 0.61103624  1  2 -0.06412938  0.17254223 -1
# 4  3 0.6110362 0.75620895  2  3 -0.06412938  0.17254223 -1
# 5  4 0.0000000 0.54577648  1  2  0.14567147  0.66730508 -1
# 6  4 0.5457765 0.74995655  2  3  0.14567147  0.66730508 -1

# The variables `Z1` and `Z2` are covariates/features to be used for tailoring treatment to individuals. 
# Estimation of an optimal individual treatment rule for prolonging the time spent in State 2 based on a linear decision function and considering the process up to time τ = 3 can be performed as follows:

fit <- itr(
  data= ret_data, 
  feat=c("age", "type", "risk"), 
  # w = c(0, 1, 0), 
  tau = 3, 
  lambda=1, 
  kernel="linear", 
  SE=TRUE
)

# print("ITR FIT")
# print(fit)

# The estimates of the coefficients of the optimal linear decision function can be obtained as follows
fit$beta_opt

# The estimated value function of the estimated individualized treatment rule can be obtained as follows
fit$V_opt

# The estimated standard error of the value function of the estimated rule can be obtained as follows
fit$se_V_opt


res = list()

# Estimation of the the value function of the latter estimated optimal treatment rule and its standard error using the function `V_d` can be performed as follows:
res$opt = V_d(data=data, w=c(0, 1, 0), tau=3, dec.fun=fit$fit, feat=c("Z1", "Z2"), SE = TRUE)

# Estimation of the the value function of the fixed rule that assigns treatment 1 to everyone, along with its standard error, can be performed as follows:
res$pos = V_d(data=data, w=c(0, 1, 0), tau=3, dec.fun=1, feat=c("Z1", "Z2"), SE = TRUE)

# Estimation of the the value function of the fixed rule that assigns treatment -1 to everyone, along with its standard error, can be performed as follows:
res$neg = V_d(data=data, w=c(0, 1, 0), tau=3, dec.fun=-1, feat=c("Z1", "Z2"), SE = TRUE)


print(res)

print("optimal")
print(fit$fit)





###########################################################################
## Estimating optimal Treatment Regime for Single Stage Trial
###########################################################################
#library(DynTxRegime) 
#library(dtrSurv)


## 1. Format data for OWL Model 
#data_ret_wl = retinopathy # [["laster", "eye", "age" ... ]]

## remove 
#data_ret_wl = retinopathy[ , !names(retinopathy) %in% c("id" )]
#data_ret_wl_y = retinopathy$futime
#data_ret_wl_tx = retinopathy$trt




## Model Single-Stage 
#dtr_survival_model = dtrSurv(
#  data = data_ret_wl,
#  txName = "trt",
#  models = Surv(futime,status) ~ laser + eye + age + type + risk
#  # tau = /C
#)


#predict(dtr_survival_model,newdata=NULL, stage = 1, findOptimal = TRUE)



## DynTxRegime 
## How to use these with Survival Time Data  ??? Can i fit a sirvival model to the buildModelObj ??? 

## # 2. Model Propensity Score 
## model_ret_main_terms =  ~ laser + eye + age + type  + risk 
## 
## # defines the model and R methods to be used to obtain parameter estimates and predictions for the propensity for treatment
## 
## # How to define this for survival outcome ? 
## moPropen <- buildModelObj(
##   model = model_ret_main_terms, 
##   solver.method = 'glm',
##   solver.args = list('family'='poisson'),
##   predict.method = 'predict.glm',
##   predict.args = list(type='response') #, propen.missing = 'smallest')
## )
## 
## # 3. outcome weighted learning 
## # tx.rules = function(a,b,c,data){
## #   as.numeric(a + b*)
## # }
## 
## owl(
##   moPropen = moPropen,
##   data = data_ret_wl,
##   reward = data_ret_wl$futime, # response vector 
##   txName = "trt",
##   regime = model_ret_main_terms, # formula obj
##   # lambdas, 
##   # cvFolds, 
##   kernel = "linear",
##   # kparam,
##   # surrogate,
##   verbose = 1
## )








