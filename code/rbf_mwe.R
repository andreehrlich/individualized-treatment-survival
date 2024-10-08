library(data.table) # data manipulation
library(radiant.data) # which.pmin()
library(survival) # methods + data sets 
library(kernlab)  # rbf kernel
library(pbapply)  # parallel
library(parallel)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
source("../models/msowl/R/beta_version/msowl_new.R")

########################################################
# work out svm rbf "prediction" (classification on new)
########################################################
msowl_rbf_s1_fit <- msowl(
  formula=ms_formula,
  id = "id",
  w = w_single_state,
  tau = tau,
  data = as.data.frame(train_data),
  kernel="rbf",
  sigma=1,
  lambda = 1,
  jackknife = FALSE,
  trim = NULL
)
msowl_rbf_fit$Value
msowl_rbf_fit$fit$alpha1
dat = test_data
feat = c("z1","z2")
dat[,c("z1","z2")]

# f(x) = sum(a_i * K(x_i, x_{new})) + b
new_dat = as.matrix(test_data[,feat])
rbf = rbfdot(sigma = msowl_rbf_fit$sigma)
K <- kernelMatrix(
  rbf,
  new_dat,
  as.matrix(msowl_rbf_fit$fit$Z)
)
f_n <- msowl_rbf_fit$fit$beta0 + K%*%msowl_rbf_fit$fit$alpha1
########################################################
