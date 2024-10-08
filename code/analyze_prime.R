# ITR analysis on real-data set: colorectal cancer panitumumab 
# PRIME Trial 
# Primary Endpoint: PFS Survival 
# Secondary Endpoint: OS Survival
# Population Strata: KRAS tumor (wild-type vs mutant)
library(data.table)
library(parallel)
library(pbapply)
library(sas7bdat) 
library(data.table)
library(survival)
library(ranger)
library(ggplot2)
library(dplyr)
library(ggfortify)
library(ggsurvfit)
library(survminer)
library(xtable)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("./msowl_modified.R")
source("./survival.regression.R")

# PRIME DATA -------
d1 = setDT(read.sas7bdat("./data/Colorec_Amgen_2006_309/PDS_DSA_20050203/adsl_pds2019.sas7bdat"))
d2 = setDT(read.sas7bdat("./data/Colorec_Amgen_2006_309/PDS_DSA_20050203/biomark_pds2019.sas7bdat"))

# exploratory data analysis 
###########################
# table(d2$BMMTR1)
# KRAS exon 2 (c12/13) 
KRAS.vec = c("Wild-type", "Mutant")
KRAS.cur = "Wild-type"
d1[d2, kras:= i.BMMTR1, on = .(SUBJID)] 
table(d1$kras)

# REMOVE PATIENTS WITHOUT KRAS INFO
d1 = d1[(kras %in% KRAS.vec), ]
d1

# SHORTEN VARIABLE NAMES
d1[kras == "Wild-type", kras:= "WT"]
d1[kras == "Mutant", kras:= "M"]
d1[TRT == "FOLFOX alone", TRT:= "F"]
d1[TRT == "Panitumumab + FOLFOX", TRT:= "P+F"]

# DEATH DAY & STATUS 
d1[, DTHDY]
d1$DTH

# Progression-Free Survival PFS -------
# PFS Day (Central, RECIST)
d1$PFSDYCR
d1$PFSCR # PD on Study

d1[, .(PFSCR, PFSDYCR, DTH,DTHDY)]
# d1[, list("ind.eq" = PFSCR == DTH, "pfs.lt.dth"= PFSDYCR <= DTHDY), by = .(PFSCR,DTH)]

# OVERALL SURVIVAL  -------
# histogram survival + censor times 
ggplot(
  d1, 
  aes(x=DTHDY, fill=as.factor(DTH))
  ) + 
  geom_histogram(
    position="dodge"
  ) + 
  facet_wrap(
    ~ kras
  ) + 
  guides(
    fill = guide_legend(title = "DTH")
  ) #+
  # ggtitle(paste0("KRAS")) +
  # xlab("Days until Death")

fname = paste0("./documents/thesis_paper/figures/prime_hist_by_kras.pdf")
# print(paste("write to file: ", fname))
# ggsave(fname, width=8, height=3,dpi=96)


ggplot(
  d1, 
  aes(x=PFSDYCR, fill=as.factor(PFSCR))
) + 
  geom_histogram(
    position="dodge"
  ) + 
  facet_wrap(
    ~ kras
  ) + 
  guides(
    fill = guide_legend(title = "PD")
  ) #+
  
# ggtitle(paste0("KRAS")) +
  # xlab("Days until Disease Progression")

fname = paste0("./documents/thesis_paper/figures/prime_hist_pfs.pdf")
# print(paste("write to file: ", fname))
# ggsave(fname, width=8, height=3,dpi=96)



## SURVIVAL PLOTS -------

# OVERALL SURVIVAL (OS)
os.fit <- survfit(Surv(DTHDY, DTH) ~ TRT + kras, data = d1)
# summary(os.fit)
ggsurvplot(
  os.fit, 
  data=d1,
  conf.int = F,       
  legend = "bottom",
  linetype=c(2,2,1,1)#"strata"
)
fname = paste0("./documents/thesis_paper/figures/prime_os_surv_no_ci.pdf")
# print(paste("write to file: ", fname))
# ggsave(fname, width=8, height=6,dpi=96)

# PROGRESSION-FREE SURVIVAL (PFS)
pfs.fit <- survfit(Surv(PFSDYCR, PFSCR) ~ TRT + kras, data = d1)
ggsurvplot(
  pfs.fit, 
  data=d1,
  conf.int = F,      
  # legend = "none",
  # linetype="strata"
  linetype=c(2,2,1,1)
)
fname = paste0("./documents/thesis_paper/figures/prime_pfs_surv_no_ci.pdf")
# print(paste("write to file: ", fname))
# ggsave(fname, width=8, height=6,dpi=96)


# BOTH ON ONE PLOT 
surv.plot.both = ggsurvplot(
  list(OS=os.fit,PFS=pfs.fit), 
  combine=TRUE,
  data=d1,
  conf.int = F ,     
  # legend = "none",
  tables.theme = theme_cleantable(),  # Clean risk table
  # palette = "lancet"
  linetype=c(2,2,1,1,2,2,1,1)
)
surv.plot.both
fname = paste0("./documents/thesis_paper/figures/prime_both_surv_no_ci.pdf")
# ggsurvplot ggsave workaround
# print(paste("write to file: ", fname))
# pdf(fname)
# print(surv.plot.both, newpage = FALSE)
# dev.off()

# TODO filter by status indicators
# median PFS by kras type
d2 = d1[PFSCR ==1 , list("Median PFS"= median(PFSDYCR/30), "Median OS"= median(DTHDY/30), "PFS Events"=.N), by=.(kras,TRT)]
d2
d2 = d1[DTH ==1 , list("Median PFS"= median(PFSDYCR/30), "Median OS"= median(DTHDY/30), "OS Events"=.N), by=.(kras,TRT)]
d2

print(
  xtable(
    d2,
     caption = c(
       "PRIME: Median Survival by Strata", 
       "PRIME: Median Survival by Strata"
     )
 ),
 row.names = F
 #, file="./documents/thesis_paper/median_pfs.tex"
)

# pan
# dt %>% filter(event.status == 1 & a == 1 ) %>% summarize(median_surv = median(t2)) / 30
# chemo
# dt %>% filter(event.status == 1 & a == -1 ) %>% summarize(median_surv = median(t2)) / 30



# ITR Estimation --------
## DATA FORMATTING -------

d1[, .N, by=.(TRT)]  # TREATMENT -1/+1 Binary values
d2 = copy(d1)
d2[TRT == "P+F", TRT := 1]
d2[TRT == "F", TRT := -1]
res_col = c("TRT")
d2[ , (res_col) := lapply(.SD, as.integer), .SDcols = res_col] # convert to integer 
d2[, AGE := scale(AGE)] # center 
d2[, B_ECOG := as.integer(B_ECOG == "Fully active")]
d2[, SEX := as.integer(SEX == "Male")]
d2[, kras := as.integer(kras == "M")]
cols.id = c("SUBJID","TRT")
col.covariates = c("AGE", "SEX","B_ECOG", "kras")
cols.keep = c(
  "SUBJID",
  "TRT",
  "AGE",
  "SEX",
  "kras",
  "B_ECOG",
  "DTH",
  "DTHDY",
  "PFSCR",
  "PFSDYCR"
)
d2 = d2[, cols.keep, with=F] 
# sd(d2$AGE)
# 2024/07/01
# TODO Multistate MSOWL
# TODO Duration in response
# - Duration of Response (DoR) 
# - Start: tumor response 
# - End: disease progression / death 
# TODO report treatment rule found by each method. 

## ITR Model Fitting ------
#  PFS Survival 

tau = 30 # try 90th-%ile ? 
month = 365.25/12
d2$PFSDYCR.m = (d2$PFSDYCR + 1) / month
d2$DTHDY.m = (d2$DTHDY + 1) / month
dt.prime = data.table(
  id = as.integer(d2$SUBJID),
  a = d2$TRT,
  # COVARIATES
  age = d2$AGE,
  sex = d2$SEX,
  kras = d2$kras,
  ecog = d2$B_ECOG,
  # times 
  t1 = 0, 
  # PFS 
  pfs.event.status = (d2$PFSCR | tau < d2$PFSDYCR.m ),
  pfs = pmin(d2$PFSDYCR.m, tau),
  pfs.s1 = 1,
  pfs.s2 = 1 + (d2$PFSCR | tau < d2$PFSDYCR.m),
  # OS 
  os.event.status = (d2$DTH | tau < d2$DTHDY.m),
  os = pmin(d2$DTHDY.m, tau), 
  os.s1 = 1,
  os.s2 = 1 + d2$DTH | tau < d2$DTHDY.m 
)

dt.prime[, .N, by=pfs.event.status ]
dt.prime[, .N, by=os.event.status ]
ggplot(dt.prime, aes(x=os, color=as.factor(os.event.status))) + geom_histogram()

pfs.formula = ms(t1 = t1, t2 = pfs,  s1 = pfs.s1, s2 = pfs.s2, a = a) ~ age + sex + ecog + kras
os.formula =  ms(t1 = t1, t2 = os,   s1 = os.s1,  s2 = os.s2,  a = a) ~ age + sex + ecog + kras
n.total = dim(dt.prime)[[1]]




# KOSOROK2015 To avoid overfitting, we employ a crossvalidated analysis. At each run, we partition the whole data set into 5 pieces, where 4 parts
# of the data are used as training data to estimate the individualized treatment rules, and the
# remaining part is the validation set for implementing the estimated rules, with empirical
# values stored for each method respectively. The cross-validated values are obtained by
# averaging the empirical values on all 5 validation subsets. The procedure is repeated 100
# times.


# PFS -------------
## MSOWL ------------
### Lambda Tuning -------------
pfs.lambda = select.lambda.kfold(
  formula=pfs.formula,
  id = "id",
  w = c(1,0),
  tau = tau,
  data = as.data.frame(dt.prime),
  # n.batch=200,
  lambdas = c(0.01, 0.05, 0.1, 0.5, 1, 5, 10, 50, 100),
  SE=T,
  debug=T
)
write.csv(
  x = pfs.lambda$details,
  "./prime_analysis/pfs_lambda_tuning_5fold.csv"
)


pfs.lambda$details = as.data.frame(pfs.lambda$details)
pfs.lambda$details$V = round(pfs.lambda$details$V,1)
pfs.lambda$details$SD = round(pfs.lambda$details$SD,1)
print(
  xtable(
    x = pfs.lambda$details,
    caption = "RCT PRIME PFS: MSOWL 5-Fold Cross Validation"
  ),
  file = "./documents/thesis_paper/prime_msowl_5fold_tuning"
)
df = data.frame(pfs.lambda$details)
ggplot(df, aes(x=log(lambdas), y=V)) + 
  geom_point() + 
  geom_errorbar(aes(log(lambdas),ymin=V-SD, ymax=V+SD)) +
  theme_minimal()

ggsave("./documents/thesis_paper/figures/prime_msowl_5fold_tuning.pdf",
       width = 4, height = 3)

### TRAIN ---------
timer = Sys.time()
print(timer)
# precalc
# msowl.reward = precalculate_reward(
#     data        = as.data.frame(dt.prime),
#     formula     = pfs.formula,
#     tau         = tau,
#     w           = c(1,0),
#     id          = "id",
#     n.batch     = 200,
#     owl_method  = "MSOWL",
#     debug       = T
# ) 
### JACKKNIFE  -------
fit.pfs = msowl(
  formula      = pfs.formula,
  id           = "id",
  w            = c(1,0),
  tau          = tau,
  data         = as.data.frame(dt.prime),
  lambda       = pfs.lambda$lambda.max*(n.total^(-1/2)),# pfs.lambda$lambda,
  jackknife    = T,
  jackknife_pg = T,
  # reward       = msowl.reward,
  debug        = T
)
fit.pfs$kernel = "linear"

# > fit.pfs$beta_opt
# (constant)           age           sex          ecog          kras 
# 2.207426e-02 -3.575219e-08  7.283798e-01  2.071013e-01 -1.064518e+00 
# > fit.pfs$beta_opt[1]
# (constant) 
# 0.02207426 
# > fit.pfs$beta_opt[2]
# age 
# -3.575219e-08 
# > fit.pfs$beta_opt[3]
# sex 
# 0.7283798 
# > fit.pfs$beta_opt[4]
# ecog 
# 0.2071013 
# > fit.pfs$beta_opt[5]
# kras 
# -1.064518 


saveRDS(fit.pfs, file="./prime_analysis/fit.pfs.RData")


res.pfs = cross.validation(
    formula = pfs.formula,
    id = "id",
    w = c(1,0),
    tau = tau,
    data = as.data.frame(dt.prime),
    lambda = pfs.lambda$lambda.1se, # 1, # pfs.lambda$lambda.1se*(n.total^(-1/2)),
    k = 5,
    debug = T,
    SE=T,
    owl_method = "MSOWL",
    B = 10
)
mean(res.pfs$V)
# [1] 11.83408
mean(res.pfs$SD)
# [1] 1.364368

# ONE-SIZE-FITS-ALL
{
  val.pfs = msowl.val(
    formula      = pfs.formula,
    id           = "id",
    w            = c(1,0),
    rule         = fit.pfs,
    tau          = tau,
    data         = as.data.frame(dt.prime),
    fixed.rules    = T,
    SE           = T,
    debug        = T
  )
  # [[1]]
  # estimate    SE     ll     ul
  # V(dn)   12.189 0.613 10.988 13.390
  # V(1)    11.424 0.479 10.485 12.363
  # V(-1)   11.083 0.370 10.358 11.808
  # [[2]]
  # estimate     ll    ul
  # V(dn) - V(1)     0.765 -0.159 1.689
  # V(dn) - V(-1)    1.106 -0.105 2.316
  # $n.failed
  # [1] 0

  # JOIN TABLES INTO ONE 
  koko = as.data.frame(val.pfs[[2]])
  koko$SE = rep(NA,2)
  koko = koko[,c(1,4,2,3)]
  dt.out = rbind(val.pfs[[1]], koko)
  dt.out = round(dt.out,1)
  print(
    xtable(
      dt.out, #as.data.frame(fit.pfs.val[[1]]),
      caption="PRIME PFS: Estimated Value (Months)"
    ),
    file="./documents/thesis_paper/prime_pfs_msowl.tex",
    hline.after = c(-1,0,3,5)
  )
}

# debug = 1
## ICO -------
debug = 1
if (debug) print("ICO")
### Lambda Tuning -------------
# ico.fit.pfs = msowl(
#   formula      = pfs.formula,
#   id           = "id",
#   w            = c(1,0),
#   tau          = tau,
#   data         = as.data.frame(dt.prime),
#   lambda       = (n.total^(-1/2)),# ico.pfs.lambda$lambda,
#   jackknife    = F,
#   jackknife_pg = T,
#   owl_method   = "ICO",
#   debug        = T
# )

ico.pfs.lambda = select.lambda.kfold(
  formula=pfs.formula,
  id = "id",
  w = c(1,0),
  tau = tau,
  data = as.data.frame(dt.prime),
  # n.batch=200,
  lambdas = c(0.01, 0.05, 0.1, 0.5, 1, 5, 10, 50, 100),
  owl_method = "ICO",
  SE=T,
  debug=T
)
write.csv(
  x = ico.pfs.lambda$details,
  "./prime_analysis/ico_pfs_lambda_tuning_5fold.csv"
)

ico.pfs.lambda$details = as.data.frame(ico.pfs.lambda$details)
ico.pfs.lambda$details$V = round(ico.pfs.lambda$details$V,1)
ico.pfs.lambda$details$SD = round(ico.pfs.lambda$details$SD,1)
print(
  xtable(
    x = ico.pfs.lambda$details,
    caption = "RCT PRIME PFS: ICO 5-Fold Cross Validation"
  ),
  file = "./documents/thesis_paper/prime_ico_5fold_tuning"
)
df = data.frame(ico.pfs.lambda$details)
ggplot(df, aes(x=log(lambdas), y=V)) + 
  geom_point() + 
  geom_errorbar(aes(log(lambdas),ymin=V-SD, ymax=V+SD)) +
  theme_minimal()

ggsave("./documents/thesis_paper/figures/prime_ico_5fold_tuning.pdf",
       width = 4, height = 3)

### TRAIN ---------
timer = Sys.time()
print(timer)
### JACKKNIFE  -------

# dt.nocens = dt.prime[pfs.event.status == 1, ]
ico.fit.pfs = msowl(
  formula      = pfs.formula,
  id           = "id",
  w            = c(1,0),
  tau          = tau,
  data         = as.data.frame(dt.prime),
  lambda       = 1.892528,
  # lambda       =  50* (nrow(dt.prime)^(-1/2)),
  # lambda       = ico.pfs.lambda$lambda.max*(nrow(dt.prime)^(-1/2)),  # ico.pfs.lambda$lambda.max
  jackknife    = F,
  jackknife_pg = T,
  owl_method   = "ICO",
  debug        = T
)
ico.fit.pfs$kernel = "linear"
saveRDS(ico.fit.pfs, file="./prime_analysis/ico.fit.pfs.RData")
# $beta_opt
# (constant)           age           sex          ecog          kras 
# -2.960033e-01 -1.339842e-07  1.260174e-01  2.970721e-02 -1.968660e-01 

# $Value
# [1] 11.03918
ico.fit.pfs

ico.res.pfs = cross.validation(
  formula = pfs.formula,
  id = "id",
  w = c(1,0),
  tau = tau,
  data = as.data.frame(dt.prime),
  lambda = 100, #*nrow(dt.prime)^(-1/2), # 1, # ico.pfs.lambda$lambda.1se*(n.total^(-1/2)),
  k = 5,
  debug = T,
  SE=T,
  owl_method = "ICO",
  B = 10
)

#       estimate  SE     LL      UL 
# V(dn)   11.803 1.547 8.772 14.835

mean(ico.res.pfs$V)
# [1] 11.03
mean(ico.res.pfs$SD)
# [1] 0.9993817

# ONE-SIZE-FITS-ALL
{
  ico.val.pfs = msowl.val(
    formula      = pfs.formula,
    id           = "id",
    w            = c(1,0),
    rule         = ico.fit.pfs,
    tau          = tau,
    data         = as.data.frame(dt.prime),
    fixed.rules    = T,
    SE           = T,
    debug        = T
  )
  # estimate    SE     ll     ul
  # V(dn)   11.083 0.426 10.249 11.918
  # V(1)    11.424 0.441 10.559 12.288
  # V(-1)   11.083 0.426 10.249 11.918
  # 
  # [[2]]
  # estimate     ll    ul
  # V(dn) - V(1)    -0.341 -1.597 0.915
  # V(dn) - V(-1)    0.000  0.000 0.000
  
  # JOIN TABLES INTO ONE 
  koko = as.data.frame(ico.val.pfs[[2]])
  koko$SE = rep(NA,2)
  koko = koko[,c(1,4,2,3)]
  dt.out = rbind(ico.val.pfs[[1]], koko)
  dt.out = round(dt.out,1)
  print(
    xtable(
      dt.out, #as.data.frame(fit.pfs.val[[1]]),
      caption="PRIME PFS: ICO Estimated Value (Months)"
    ),
    file="./documents/thesis_paper/prime_pfs_ico.tex",
    hline.after = c(-1,0,3,5)
  )
}

## COX -------
if (debug) print("COX")
# pfs.cox.formula = as.formula("Surv(time=pfs, pfs.event.status, type='right') ~ age + sex + ecog + kras + a + a:age + a:sex + a:ecog + a:kras")
cox.res.pfs = cross.validation(
  formula = pfs.formula,
  data = as.data.frame(dt.prime),
  k = 5,
  debug = T,
  owl_method = "COX",
  B = 10
)


# cox.itr <- function(train_set, test_set)
# FIT MODEL
ids = sort(unique(dt.prime$id))
cox.V = rep(NA,length(ids))
cox.D = rep(NA,length(ids))
for (i in seq_along(ids)){
  
  if (debug) {print(paste0("ID: ", i, "/", length(ids)) )}
  train_set = as.data.frame(dt.prime[which(id != ids[i]), ])
  test_set  = as.data.frame(dt.prime[which(id == ids[i]), ])
  
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
  pred_pos = summary(res)$table["median"]
  # -1
  data.copy$a = -1  
  res = survfit(cox.train, newdata= data.copy)
  pred_neg = summary(res)$table["median"]
  if (any(is.na(pred_pos)) | any(is.na(pred_pos))) {
    print("NA in PREDICTION! ")
    print(head(data))
    print("NA in PREDICTION! ")
  }
  better_treatment = ( which.pmax(pred_neg, pred_pos) - 1 ) * 2 - 1
  cox.test.predictions = data.frame(test_set$id, A_hat=better_treatment)
  cox.D[i] = cox.test.predictions$A_hat[1]
}
# Calculate Value
creward = calculate_reward(
  pfs.formula, 
  id="id", 
  w = c(1,0), 
  tau = tau, 
  data = dt.prime,
  owl_method = "MSOWL"
)

Wt = creward[[2]]$Wt
A = creward[[2]]$A
vecVal=Wt*(A*cox.D>=0)
V_jack <- mean(vecVal,na.rm=TRUE)
# 12.10195
SD_jack=sd(vecVal,na.rm=TRUE)
# 18.15387

########## 
train_set = as.data.frame(dt.prime) #[which(id != ids[i]), ])
test_set  = as.data.frame(dt.prime) #[which(id == ids[i]), ])

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
data.copy$a = as.factor(1)
res = survfit(cox.train, newdata= data.copy)
pred_pos = summary(res)$table[, "median"]
# -1
data.copy$a = as.factor(-1)
res = survfit(cox.train, newdata= data.copy)
pred_neg = summary(res)$table[, "median"]
if (any(is.na(pred_pos)) | any(is.na(pred_pos))) {
  print("NA in PREDICTION! ")
  print(head(data))
  print("NA in PREDICTION! ")
}
better_treatment = ( which.pmax(pred_neg, pred_pos) - 1 ) * 2 - 1
cox.test.predictions = data.frame(test_set$id, A_hat=better_treatment)

val <- msowl.val(
  formula=pfs.formula,
  id="id",
  w = c(1,0),
  tau=tau,
  rule = cox.test.predictions,
  fixed.rules = T,
  SE=T,
  debug = T,
  data = test_set
)
if (debug) print(val)
# V(dn)     V(1)    V(-1) 
# 12.41048 11.42381 11.08316 


# Call:
#   coxph(formula = pfs.cox.formula, data = train_set)
# n= 866, number of events= 771 
#                             coef   exp(coef)  se(coef)  z    Pr(>|z|)   
# age                      -0.02778   0.97260  0.05007 -0.555  0.57903   
# factor(sex)1              0.08314   1.08670  0.10855  0.766  0.44371   
# factor(ecog)1            -0.26813   0.76481  0.10257 -2.614  0.00894 **
# factor(kras)1             0.06601   1.06823  0.10383  0.636  0.52494   
# factor(a)1                0.10417   1.10979  0.16428  0.634  0.52601   
# age:factor(a)1            0.06572   1.06792  0.07588  0.866  0.38649   
# factor(sex)1:factor(a)1  -0.20436   0.81517  0.15481 -1.320  0.18681   
# factor(ecog)1:factor(a)1 -0.23189   0.79304  0.14631 -1.585  0.11300   
# factor(kras)1:factor(a)1  0.35083   1.42024  0.14735  2.381  0.01727 * 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
#                            exp(coef) exp(-coef) lower .95 upper .95
# age                         0.9726     1.0282    0.8817    1.0729
# factor(sex)1                1.0867     0.9202    0.8784    1.3443
# factor(ecog)1               0.7648     1.3075    0.6255    0.9351
# factor(kras)1               1.0682     0.9361    0.8715    1.3093
# factor(a)1                  1.1098     0.9011    0.8043    1.5314
# age:factor(a)1              1.0679     0.9364    0.9203    1.2392
# factor(sex)1:factor(a)1     0.8152     1.2267    0.6018    1.1041
# factor(ecog)1:factor(a)1    0.7930     1.2610    0.5953    1.0564
# factor(kras)1:factor(a)1    1.4202     0.7041    1.0640    1.8958
# 
# Concordance= 0.588  (se = 0.011 )
# Likelihood ratio test= 49.75  on 9 df,   p=1e-07
# Wald test            = 50.6  on 9 df,   p=8e-08
# Score (logrank) test = 51.07  on 9 df,   p=7e-08


summary(cox.train)



##############

# test_result = c(
#   n     = n.train,  
#   scen  = scenario,
#   cens_param = censor,
#   pct_cens =round(pct_censored_train,2),
#   b     = b,
#   dur   = round(Sys.time() - round_time,2),
#   PLUS  = unname(msowl_test["V(1)"]),
#   MINUS = unname(msowl_test["V(-1)"]),
#   MSOWL = unname(msowl_test["V(dn)"]),
#   ICO   = unname(ico_test),
#   COX   = unname(cox.test.val["V(dn)"])
# )


