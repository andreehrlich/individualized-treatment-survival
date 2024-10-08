# ZhaoKosorok2015 ITR Simulation Experiment Replication
# MSc Statistics Thesis: Individualized Treatment Regimes (ITR)  2024
# Andre Ehrlich
# Supervisors: Demiris, Bakoyannis
library(pbapply)
pbo = pboptions(type="txt") # timer
library(parallel)
library(radiant.data) # which.pmin()
library(ggh4x)
library(ggplot2)
# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

source("./msowl_modified.R")
# source("./msowl_utility.R")
source("./sim_surv_regr.R")
source("./sim_kos.R")

# MAIN SIM LOOP -----
doOne.div = function(my_input) {
  round_time = Sys.time()
  scenario   = my_input[[1]]
  censor     = my_input[[2]]
  n.train    = my_input[[3]]
  n.test     = my_input[[4]]
  n.batch    = my_input[[5]]
  tau        = my_input[[6]]
  b          = my_input[[7]]
  debug      = my_input[[8]]
  div.method = my_input[[9]]
  lambda     = my_input[[10]]
  ico.lambda = my_input[[11]]
  # data.set   = my_input[[12]]

  if (debug) {print("SIMULATION LOOP"); print(my_input)}
  if (scenario <= 4) {
    formula = ms(t1 = t1, t2 = t2, s1 = s1, s2 = s2, a = A) ~ z1 + z2 + z3
    w = c(1,0)
    # TRAIN  -----------------------------------------------------------
    data.train = sim.kos(
      n = n.train,
      scenario = scenario,
      tau = tau,
      censor.rate = censor
    )
    # TEST  -----------------------------------------------------------
    data.test = sim.kos(
      n = n.test,
      scenario = scenario,
      tau = tau,
      censor.rate = FALSE
    )
  } else if (scenario == 5) {
    # There may be a problem with my data generation function (DGF)
    # Test with Bakoyannis's provided example data.
    # A few considerations:
    #  - it is multi-state
    #     --> How to mimic simple survival ? 
    #  - N = 200, with 291 total rows for multi-state data. 
    #     --> How to synthesise more examples from this data?  SMOTE ? 
    formula = ms(t1 = t1, t2 = t2, s1 = s1, s2 = s2, a = A) ~ z1 + z2 
    w = c(1,0,0)
    {
      data.example = setDT(read.csv("./models/msowl/data/example_data.csv"))
      # data.synth = c()
      # for (i in 1:100){
      #   data.example.2 = copy(data.example)
      #   data.example.2[, 
      #     `:=`(
      #       z2 = z2 + rnorm(nrow(data.example.2),0,.05), 
      #       id = id + max(data.example$id)*i
      #     )
      #   ]
      #   data.synth = rbind(data.synth, data.example.2)
      # }
      # data.synth = setDT(data.synth)
      data.example[, status := as.numeric(s1 != s2 & t2 <= tau) ]
      data.train = data.example[1:n.train,]
      data.test = data.example[n.train:(n.test + n.train-1),]
    }
  }
  
  # print(data.test)
  # CENSORING PERCENTAGE ----------------
  train.pct.censored = mean(!data.train$status)
  test.pct.cens = mean(!data.test$status)
  if (debug) {
    print(paste("train cens pct: ", train.pct.censored))
    print(paste("test cens pct: ",test.pct.cens))
  }
  
  # PRECALC TEST ----------------------------------------------------
  time.pre.calc = Sys.time()
  data.test.complete = FALSE
  if (div.method == 2) {
    data.test.complete = precalculate_reward(
      data=data.test,
      id = "id",
      w = w,
      formula=formula,
      tau=tau,
      n.batch = n.batch
    )
  }
  time.pre.calc = Sys.time() - time.pre.calc
  if (debug) print(paste("PRE CALC: ", time.pre.calc))
  
  
  ## MSOWL-----------------------------------------------------------
  ### MSOWL TRAIN---------
  if (debug) print("MSOWL")
  msowl_train <- msowl(
    formula = formula,
    id = "id",
    w = w,
    tau = tau,
    data = as.data.frame(data.train),
    lambda = dim(data.train)[[1]]^(-1/2) * lambda,
    jackknife = FALSE,
    trim = NULL,
    reward = NULL,
    debug=F
  )
  if (debug) {print("msowl_train"); print(msowl_train)}
  if(unlist(gregexpr('error', class(msowl_train)))[1] !=-1){ 
    # ERROR
    itr.msowl.val = list('V(dn)'=NA, "V(1)"=NA, "V(-1)"=NA)
    time.msowl.val = NA
  } else{
    ### MSOWL TEST ------------
    if (debug) print(msowl_train$Value)
    time.msowl.val = Sys.time()
    if (div.method == 1) {
      batches = split(
        as.data.frame(data.test),
        sample.int(
          n=n.test/n.batch,
          size=nrow(data.test),
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
        w = w,
        tau = tau,
        fixed.rules=T,
        reward=data.test.complete
      )
    } else if (div.method == 3){
      itr.msowl.val = msowl.val(
        data = as.data.frame(data.test),
        rule = msowl_train,
        formula = formula,
        id = "id",
        w = w,
        tau = tau,
        fixed.rules=T,
        reward=NULL
      )
    }
    time.msowl.val = difftime(Sys.time(), time.msowl.val, units = "secs")
    if (debug) print(paste(c("MSOWL", itr.msowl.val, time.msowl.val)))
  }
  ## AGG RESULTS ----------
  test_result = c(
    n.train       = n.train,  
    n.test        = n.test,
    scen          = scenario,
    cens_param    = censor,
    pct.cens      = round(train.pct.censored,2),
    test.pct.cens = test.pct.cens,
    b             = b,
    dur           = difftime(Sys.time(), round_time, units = "secs"), 
    COMP_METHOD   = div.method,
    lambda        = lambda,
    PLUS          = unname(itr.msowl.val["V(1)"]),
    MINUS         = unname(itr.msowl.val["V(-1)"]),
    MSOWL         = unname(itr.msowl.val["V(dn)"]),
    MSOWL.dur     = time.msowl.val
  )
  if (debug) print(test_result)
  return(test_result)
}



# SIM: VALUE FUNC SAMPLE SIZE --> PROC. TIME -----------
print("VALUE FUNC: SAMPLE SIZE --> PROC. TIME")
div.method.dict = list(
  "1" = "DIV.SIMPLE",
  "2" = "PRECALC",
  "3" = "NO.DIV"
)

cl = makeCluster(
  detectCores() - 1,
  outfile="./sim_log.txt"
)

clusterExport(
  cl,
  varlist = c(
    "copy",
    "sim.kos",
    "inv.trans.sample.cox",
    "H.inv.kos",
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
    "select.lambda.kfold",
    "pbsapply"
    # "select.lambda.new"
  )
)



compare.div.methods = function(
  n.test.vec = c(500,1000,1500,2000),
  B = 20
) {
  
  sim.timer = Sys.time()
  agg.results = data.frame()
  for (n.test in n.test.vec) {
    for (censor in c(F,-1)) {
      for (div.method in c(2,3)) {
        print(Sys.time())
        print(paste0(n.test,censor,div.method))
        my_inputs = data.frame(
          scenario = rep(1,B),
          censor = rep(censor,B),
          n.train =  rep(400,B),
          n.test =  rep(n.test,B),
          n.batch =  rep(200,B),
          tau =  rep(tau.vec[1],B),
          b =  1:B,
          debug = rep(F,B),
          div.method =  rep(div.method,B),
          lambda =  rep(0.01,B),
          ico.lambda = rep(0.01 ,B)
        )
        
        ### PBAPPLY -----
        print("INPUTS")
        print(my_inputs)
        print(dim(my_inputs))
        print(paste("parrallel start time:", Sys.time()))
        cres = pbapply(
          X = my_inputs, 
          FUN=doOne.div, 
          MARGIN=1,
          cl= cl
        )
        
        print(t(cres))
        ### WRITE TO FILE -------
        c_fname = paste0("./sim_results/div_compare_n.test=",n.test, "_censor=",censor,"div.method=",div.method,".csv")
        print(paste("write to file.", c_fname))
        write.table(
          t(cres),
          c_fname,
          append = FALSE,
          row.names=FALSE,
          col.names=TRUE,
          sep=","
        )
        
        agg.results = rbind(agg.results, t(cres))
      }
    }
  }
  stopCluster(cl)
  print(agg.results)
  print(paste0("Duration: ", Sys.time() - sim.timer))
  return(agg.results)
}

res = compare.div.methods()
res = read.csv("./sim_results/div_compare/combined.csv")
res = setDT(res)
hist(res$dur)
res[
  , 
  list(msowl = mean(MSOWL),plus= mean(PLUS), minus= mean(MINUS),B=max(b)), 
  by=.(COMP_METHOD,cens_param, n.test)
]
id_vars = c("n.train", "n.test", "scen", "cens_param", "pct.cens","b","dur", "lambda", "COMP_METHOD")
mvars = c("PLUS","MINUS","MSOWL")
dt_long <- melt(
  res, 
  measure.vars = mvars,
  id.vars= id_vars
)
dt_long$cens.title = "Cens Param"
dt_long$model.title = "Model"
dt_long$
# dt_long[scen %in% c(1,2), surv.distr := "COX"] 
# dt_long[scen %in% c(3,4), surv.distr := "AFT"] 
# dt_long[scen %in% c(5,6,7,8), surv.distr := "EXP"] 
# dt_long[scen %in% c(9), surv.distr := "COX-NL"] 
# dt_long$surv.distr = factor(dt_long$surv.distr, levels = c("COX","AFT","EXP","COX-NL"))
ggplot(
  dt_long, 
  aes(x=factor(n.test), y=log(value), color=as.factor(COMP_METHOD))
) + 
  geom_boxplot() +
  facet_nested(
    cens.title + cens_param ~ model.title + variable
    # scales = "free"
  ) + 
  # theme_minimal() +
  theme(
    # axis.title.x=element_blank(),
    # axis.ticks.x=element_blank(),
    legend.position="bottom"
  ) + ggtitle(
    paste0(
      "B=", max(dt_long$b)
      # "; n.train=",unique(dt.kos$n),
      # "; ico.na.rate=", ico.na.rate,
      # "; msowl.na.rate=",msowl.na.rate,
      # "; Censoring Rate: ", round(mean(dt_long$pct.cens),2),
      # "; lambda=", unique(dt.kos$lambda)
    )
  )
# SAVE TO FOLDER FOR LATEX REPORT
fname = paste0("./documents/thesis_paper/figures/div.conq.value.compare.pdf")
print(paste("saving plot to file: ", fname))
ggsave(fname, width=8, height = 6)


title = "msowl_val_time"
fname.csv = paste0("./sim_results/", title, ".csv")
write.csv(x = res, file = fname.csv)
print(paste0("write csv to file: ", fname.csv))

# ggplot( lala$n.test, lala$MSOWL.dur)
ggplot(
  res, 
  aes(x=n.test,y=MSOWL.dur)
)  + geom_point() +
  + geom_line(aes(x=x,y=y_hat)) +
  theme(
    legend.position = "none"
  ) + labs(
    y  = "Duration (Seconds)",
    x = "Sample Size (N)"
  )


fname = paste0("./documents/thesis_paper/figures/",title,".pdf")
print(fname)
ggsave(
  fname,
  width=8,
  height = 6
)

# CACHED RESULTS
res = data.frame(
  x = c(500,1000,1500,2000,5000,10000),
  no.div = c(   1.9,  9.6, 25, 52, 612, 5242),
  div.conq = c( 1.9,    4,  6, 10,  20, 27)
)

# no.div cubic fit
fit = lm(no.div ~ I(x^3), data=res)
summary(fit)
linedat = data.frame(x=1:10000)
linedat$no.div.cubic = predict(fit,linedat)
tail(linedat)

# div.conq linear fit 

fit = lm(div.conq ~ x, data=res)
summary(fit)
linedat$div.conq.linear = predict(fit,linedat)
tail(linedat)


ggplot()  +
  geom_point(data = res,     aes(x = x, y = no.div,     color = "red")) +
  geom_point( data = res, aes(x = x, y = div.conq, color = "green")) +
  geom_line( data = linedat, aes(x = x, y = no.div.cubic,   color = "red"  ), linetype="dotted") +
  geom_line( data = linedat, aes(x = x, y = div.conq.linear, color = "green"), linetype="dotted") +
  scale_colour_manual(
    name   = 'Method',
    values = c(
      'red' = 'red', 
      "green" = "green",
      'red' = 'red',
      "green" = "green"
    ),
    labels = c(
      'Div&Conq',
      'No Div',
      'Cubic Fit',
      "Linear Fit"
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title=element_blank()
  ) +
  labs(
    y     = "Duration (Seconds)", 
    x     = "Sample Size (N)"
    # title = "Sample Size x Processing Time"
  )

title = "approx_msowl_time"
fname = paste0("./documents/thesis_paper/figures/",title,".pdf")
print(fname)
ggsave(
  fname,
  width=8,
  height = 6
)

# # RUN SIMULATION ---------
# run.sim.kos(
#   keyword = "dq1",
#   B = 100,
#   n.vec = c(200),
#   n.test = 400,
#   censor.vec = c(F,-1),
#   scenario.vec = c(2)
# )
#######################
# # MC SIM  -----
# run.sim.kos <- function(
    #     B = 1000,
#     n.test = 10000,
#     n.batch = 200,
#     n.vec = c(100, 200, 400),
#     scenario.vec = c(1, 2, 3, 4),
#     censor.vec = c(-1, -0.5, -0.2),
#     lambda.vec = NULL,
#     keyword = "KOS_SIM",
#     kfolds = 10
#     # num_cores = detectCores() - 1
# ) {
#   print("KOS SIMULATION")
#   print("B")
#   print(B)
#   print("n.test")
#   print(n.test)
#   print("n.batch")
#   print(n.batch)
#   print("scenario.vec")
#   print(scenario.vec)
#   print("censor.vec")
#   print(censor.vec)
#   print("lambda.vec")
#   print(lambda.vec)
#   print("n.vec")
#   print(n.vec)
#   print("keyword")
#   print(keyword)
#   
#   print("kfolds")
#   print(kfolds)
#   
#   ## PARALLEL  -----
#   # Parallel Number of Cores preferences between local mac VS remote windows machines
#   switch(
#     Sys.info()[['sysname']],
#     Windows= {
#       # remote compute 
#       print("I'm a Windows PC. Use all cores")
#       num_cores = detectCores()
#       print(paste("num_cores", num_cores))
#     },
#     Darwin = {
#       # personal computer
#       print("I'm a Mac. Use all cores but one.")
#       # source("./models/msowl/R/beta_version/msowl_new.R")
#       num_cores = detectCores()-1
#       print(paste("num_cores", num_cores))
#     }
#   )
#   # WINDOWS MACHINE CONFIGURATION
#   # cl - A cluster object created by makeCluster, 
#   #    - OR an integer to indicate number of child-processes 
#   # integer values are ignored on Windows) for parallel evaluations 
#   cl = makeCluster(
#     num_cores,
#     outfile=""
#   )
#   
#   clusterExport(
#     cl,
#     varlist = c(
#       "copy",
#       "sim.kos",
#       "inv.trans.sample.cox",
#       "H.inv.kos",
#       "precalculate_reward",
#       "calculate_reward",
#       "coxph",
#       "Surv",
#       "survfit",
#       "ms",
#       "zeta_t",
#       "basehaz",
#       "bind_rows",
#       "msowl",
#       "msowl.val",
#       "wsvm_solve",
#       "ipop",
#       "primal",
#       "coxph_itr_predict",
#       "which.pmax",
#       "select.lambda.kfold"
#       # "select.lambda.new"
#     )
#   )
#   
#   # save to folder.
#   foldername = paste0("./sim_results/", keyword)
#   dir.create(foldername, showWarnings = FALSE)
#   saved_files = c()
#   for (n.train in n.vec){
#     for (censor in censor.vec) {
#       for(scenario in scenario.vec){
#         ## FNAME ------
#         c_fname = paste0(foldername,"/c",censor,"_s",scenario,"_n", n.train,".csv")
#         print(c_fname)
#         if (file.exists(c_fname)){
#           print("this file already exists. skipping.")
#           next
#         }
#         saved_files = c(saved_files,c_fname)
# 
#         my_inputs  = data.frame(
#           scenario = rep(scenario,B),
#           censor   = rep(censor,B),
#           n.train  = rep(n.train,B),
#           n.test   = rep(n.test,B),
#           n.batch  = rep(n.batch,B),
#           tau      = rep(tau.vec[scenario],B),
#           b        = 1:B,
#           debug    = rep(F,B),
#           div.method = rep(2,B),
#           lambda   = rep(0.1,B),
#           ico.lambda = rep(0.1, B)
#           # formula = rep(my_formula,B)
#         )
# 
#         ## PBAPPLY -----
#         print("INPUTS")
#         print(my_inputs)
#         print(dim(my_inputs))
#         print(paste("parrallel start time:", Sys.time()))
#         cres = pbapply(
#           X = my_inputs, 
#           FUN=doOne.div, 
#           MARGIN=1,
#           cl= cl
#         )
#         ## FOUT -------
#         print("write to file.")
#         write.table(
#           t(cres),
#           c_fname,
#           append = FALSE,
#           row.names=FALSE,
#           col.names=TRUE,
#           sep=","
#         )
#       } # scen
#     } # cens
#   } # n.train
#   stopCluster(cl)
#   # TODO combine saved_files into single csv 
# }

#################

