# svm.lambda.tuning.R
source("./msowl_modified.R")
source("./sim_kos.R")
library(pbapply)
library(parallel)
library(data.table)
library(xtable)
library(ggh4x)
library(ggplot2)

# SVM Learning Rate Parameter --------------------------------------------------
# The learning rate parameter can be tuned according to grid search using the
# k-fold cross validation, with jackknife estimator of predictive accuracy.
# Although this is time-consuming depending on the number of possible parameter
# values chosen for comparison. $\lambda$ ought to vanish as $n \rightarrow
# \infty$, so we can safely choose $\lambda = n^{-1/2}$ which satisfies
# regularity conditions. 


# Possible Experiment:
# Do we have a notion of how this parameter really affects the estimator? 
# experiment with tuning of learning-rate/convergance parameter, and see the 
# resultant properties of the estimator with simulations. 
# $n \in 100,200,300,400,500$, $\lambda \in [10^-3, 10^3]$ in 12 steps. 
# 
# In a simulation experiment where there are many sample sizes, censoring rates,
# and linear predictors specified in decision-making statistical model fitting.
# I want to tune the sim underlying each model before the B=1000 repetitions. Why
# do 5-fold versus jackknife versus a purely generative Monte Carlo approach ?
#   If I do 5-fold CV, do not I optimise lambda for a 1/k sized sample of that
# which will be simulated and evaluated?
#   
#   An ill-specified learning rate parameter $\lambda_n$ can lead to suboptimal or
# even non-sensical estimates, because the loss function is over or under
# penalised, the model is not fitted properly, and estimates may be inconsistent.
# The variability of the estimates can grow significantly. 
--------------------------------------------------------------------------------

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


# PLOT SIM RESULTS -----------------------------
# LAMBDA COMPARISON
{
  dt.lp1 = load.files.from.dir(keyword="fixed")
  dt.lp1$lambda = 1
  dt.lp1$pct.cens = dt.lp1$pct_cens
  
  dt.lambda = rbind(
    dt.lp1, 
    load.files.from.dir(keyword="june1"),
    load.files.from.dir(keyword="june1-lambda0.001"),
    load.files.from.dir(keyword="lambda=0.0001"),
    fill=T
  )
  
  plot.lambda(dt.lambda)
}

plot.lambda <- function(
    dt.kos,
    cens_param=-1,
    n.train=100
){
  
  dt.kos = dt.kos[cens_param == cens_param & n == n.train, ]
  
  # AVG PCT CENSORED
  dt.kos[, cens.mean := round(mean(pct.cens),2)  , by = .(cens_param)]
  
  # NA COUNT
  ico.na.rate   = paste0(round(mean(is.na(dt.kos$ICO)), 2)*100, "%")
  ico.na.rate
  msowl.na.rate = paste0(round(mean(is.na(dt.kos$MSOWL)), 2)*100, "%")
  msowl.na.rate
  
  # NA COUNT for ICO Model (by N, Scenario)
  icona = dt.kos[
    , 
    list(ico_na_rate= mean(is.na(ICO))), 
    by=.(n, scen, cens.mean)
  ]
  
  # ICO.NA x N.sample.size
  ggplot(icona, 
         # aes(n, ico_na_rate)
         aes(x=factor(cens.mean), y=ico_na_rate)
  ) + 
    geom_point(aes(colour=n), alpha=0.7) + 
    facet_wrap(~scen) +
    theme(
      legend.position="bottom"
    ) + 
    xlab("Avg. Censoring Rate") + 
    ylab("Proportion of Numerical Errors")
  # ggsave("./documents/thesis_paper/figures/ico_na_by_cens.pdf")
  
  # ICO.NA x Cens.params
  ggplot(
    icona, 
    aes(x=factor(n), y=ico_na_rate, color=factor(cens.mean))
  ) + 
    geom_point() +#es(colour=factor(scen)), alpha=0.7) + 
    facet_wrap(~scen) +
    theme(
      legend.position="bottom"
    ) + 
    ylab("Proportion of Numerical Errors") + 
    xlab("N (Sample Size)")
  # ggsave("./documents/thesis_paper/figures/ico_na_by_sample_size.pdf")
  
  # MELT
  id_vars = c("n", "scen", "cens_param", "pct.cens", "cens.mean","b","dur", "lambda")
  mvars = c("PLUS","MINUS","MSOWL", "ICO", "COX")
  dt_long <- melt(
    dt.kos, 
    measure.vars = mvars,
    id.vars= id_vars
  )
  # TITLES 
  dt_long$cens.title = "Avg Pct Censored"
  dt_long[scen %in% c(1,2), surv.distr := "COX"] 
  dt_long[scen %in% c(3,4), surv.distr := "AFT"] 
  dt_long[scen %in% c(3,4), surv.distr := "EXP"] 
  dt_long[scen %in% c(9), surv.distr := "COX-NL"] 
  # dt_long$surv.distr = factor(dt_long$surv.distr, levels = c("COX","AFT"))
  dt_long$lambda.title = "Lambda"
  # dt_long$n.label = "Sample Size"
  # BIG GRID BOXPLOT
  ggplot(
    dt_long, 
    aes(x=factor(n), y=log(value), color=variable)
  ) +  
    geom_boxplot(
    ) + 
    facet_nested(
      surv.distr + scen  ~ lambda.title + lambda,
      scales = "free"
    ) +
    theme(
      axis.title.x=element_blank(),
      axis.ticks.x=element_blank(),
      legend.position="bottom"
    ) + ggtitle(
      paste0(
        "B=", max(dt.kos$b), 
        "; n.train=",unique(dt.kos$n),
        "; ico.na.rate=", ico.na.rate,
        "; msowl.na.rate=",msowl.na.rate,
        "; Censoring Rate: ", round(mean(dt_long$pct.cens),2),
        "; lambda=", unique(dt.kos$lambda)
      )
    )
  # SAVE TO FOLDER FOR LATEX REPORT
  fname = paste0("./documents/thesis_paper/figures/kos_lam_n100_many_cens_", keyword, ".pdf")
  print(paste("saving plot to file: ", fname))
  # ggsave(fname, width=14, height = 10)
}
#######



# KFOLD SIM -----------
assess.k.fold.msowl = function(my_input) {
  round_time = Sys.time()
  scenario = my_input[[1]]
  censor   = my_input[[2]]
  n.train  = my_input[[3]]
  n.test   = my_input[[4]]
  n.batch  = my_input[[5]]
  tau      = my_input[[6]]
  b        = my_input[[7]]
  debug    = my_input[[8]]
  comp_method = my_input[[9]]
  my_formula = ms(t1 = t1, t2 = t2, s1 = s1, s2 = s2, a = A) ~ z1 + z2 + z3
  # print(my_input)
  if (debug) {print("SIMULATION LOOP"); print(my_input)}
  ## TRAIN -----
  data.train = sim.kos(
    n = n.train,
    scenario = scenario,
    tau = tau,
    censor.rate = censor
  )
  train.pct.censored = mean(!data.train$status)
  if (debug) print(paste("train cens pct: ", train.pct.censored))
  ## TEST -------------
  data.test = sim.kos(
    n = n.test,
    scenario = scenario,
    tau = tau,
    censor.rate = FALSE
  )
  test.pct.cens = mean(!data.test$status)
  if (debug) print(paste("test cens pct: ",test.pct.cens))

  # PRECOMPUTE REWARD FOR "LARGE" DATA SETS  -----------------------------------
  data.train.complete = NULL
  data.train.complete.ico = NULL
  data.test.complete = FALSE
  batches = split(
    as.data.frame(data.test),
    sample.int(
      n=n.test/n.batch,
      size=nrow(data.test),
      replace=TRUE
    )
  )
  time.pre.calc = NA
  count.pre.calc = 0
  if (comp_method == 2) {
    time.pre.calc = Sys.time()
    count.pre.calc = 0
    if(n.train >= 3*n.batch){
      # if (debug)
      count.pre.calc = count.pre.calc + 1
      data.train.complete = precalculate_reward(
        data=data.train,
        n=n.train,
        my_formula=my_formula,
        tau=tau,
        n.batch = n.batch
      )
      # print("data.train.complete")
      # print(data.train.complete)
    }
    if(n.train >= 3*n.batch){
      count.pre.calc = count.pre.calc + 1
      data.train.complete.ico = precalculate_reward(
        data=data.train,
        n=n.train,
        my_formula=my_formula,
        tau=tau,
        n.batch = n.batch,
        owl_method="ICO"
      )
    }
    # TEST
    if(n.test >= 3*n.batch){
      count.pre.calc = count.pre.calc + 1
      data.test.complete = precalculate_reward(
        data=data.test,
        n=n.test,
        my_formula=my_formula,
        tau=tau,
        n.batch = n.batch
      )
    }
    time.pre.calc = Sys.time() - time.pre.calc
    if (debug) print(paste("PRE CALC: ", time.pre.calc, count.pre.calc))
  }
  ## MODEL: MSOWL -------------------------------------------------------------
  if (debug) print("MSOWL")
  # svm.lambda = n.train^(-1/2) # regularity conditions
  
  
  lambda.tuning = select.lambda.kfold(
    formula = my_formula,
    id = "id",
    w = c(1,0),
    tau = tau,
    data = as.data.frame(data.train),
    lambdas = c(0.01, 0.1), # , 0.5, 1, 5, 10, 50, 100),
    k = 5
  )
  lambda.max = lambda.tuning$lambda
  
  msowl_train <- msowl(
    formula=my_formula,
    id = "id",
    w = c(1,0),
    tau = tau,
    data = as.data.frame(data.train),
    # kernel="linear",
    lambda = dim(data.train)[[1]]^(-1/2) * lambda.max,
    jackknife = FALSE,
    trim = NULL,
    # reward = data.train.complete,
    debug=F
  )
  if (debug) print(msowl_train$Value)
  ### TEST ---------
  if(unlist(gregexpr('error', class(msowl_train)))[1] !=-1){ # ERROR
    itr.msowl.val = list('V(dn)'=NA, "V(1)"=NA, "V(-1)"=NA)
    time.msowl.val = NA
  } else{
    time.msowl.val = Sys.time()
    if (comp_method == 1) {
      res = pbapply::pbsapply(
        batches,
        msowl.val,
        formula=my_formula,
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
    } else if (comp_method == 2) {
      itr.msowl.val = msowl.val(
        data = as.data.frame(data.test),
        rule = msowl_train,
        formula = my_formula,
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
  
  ## agg res -------------------------------------------------------------
  test_result = c(
    n     = n.train,  
    scen  = scenario,
    cens_param = censor,
    pct.cens =round(train.pct.censored,2),
    b     = b,
    dur   = round(Sys.time() - round_time,2),
    COMP_METHOD = comp_method,
    PLUS  = unname(itr.msowl.val["V(1)"]),
    MINUS = unname(itr.msowl.val["V(-1)"]),
    MSOWL = unname(itr.msowl.val["V(dn)"]),
    MSOWL.dur = time.msowl.val
  )
  if (debug) print(test_result)
  return(test_result)
}

# DEBUG -------
my_input= list(
  scenario    = 1,
  censor      = -1.4,
  n.train     = 400,
  n.test      = 10000,
  n.batch     = 200,
  tau         = tau.vec[4],
  b           = 1,
  debug       = T,
  comp_method = 2
)

# run once ------- 
res_debug = assess.k.fold.msowl(
  my_input = my_input
)
print(res_debug)

# RUN SIM -------------------------------------------------------------
run_simulation <- function(
    scenario = 1,
    censor = -1.4,
    n.train = 400,
    n.test = 10000,
    n.batch = 400,
    tau = tau.vec[4],
    debug = T,
    comp_methods = c(1,2),
    B = 2
) {
  
  cl = makeCluster(
    num_cores,
    outfile=""
  )
  clusterExport(
    cl,
    varlist = c(
      "compare.value.computation.approaches",
      "itr.cox.divide.conquer.value",
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
      "select.lambda"
      # "select.lambda.new"
    )
  )
  
  # save to folder.
  foldername = "./sim_results/val_comp"
  dir.create(foldername, showWarnings = FALSE)
  # saved_files = c()
  
  # parellel in sequential chunks to save results intermittently :) 
  c_fname = paste0(foldername,"/c",my_input$censor,"_s",my_input$scenario,"_n", my_input$n.train, ".csv")
  print(c_fname)
  
  # if (file.exists(c_fname)){
  #   print("this file already exists. skipping.")
  #   next
  # }
  # saved_files = c(saved_files,c_fname)
  
  ## inputs --------
  res.agg = list()
  for (comp_method in comp_methods) {
    print(comp_method)
    my_inputs = data.frame(
      scenario = rep(scenario,B),
      censor   = rep(censor,B),
      n.train  = rep(n.train,B),
      n.test.  = rep(n.test,B),
      n.batch  = rep(n.batch,B),
      tau      = rep(tau,B),
      b        = 1:B,
      debug    = rep(F,B),
      comp_method = rep(comp_method, B)
    )
    
    cres = pbapply(
      X = my_inputs, 
      FUN = compare.value.computation.approaches, 
      MARGIN = 1,
      cl = cl
    )
    print("cres")
    print(cres)
    print("do.call(rbind,cres)")
    cres = do.call(rbind, cres)
    print(cres)
    
    ## write to file -------
    print(paste("write to file:", c_fname))
    write.table(
      cres,
      c_fname,
      append = comp_method != 1, # append after first iteration
      row.names=FALSE,
      col.names=comp_method == 1,
      sep=","
    )
    
    # save to list in-memory
    res.agg[[comp_method]] = cres
  }
  
  # combine all tables
  res = do.call(rbind, res.agg)
  # res <- data.table::rbindlist(res.agg)
  
  stopCluster(cl)
  
  return(res)
}


# RUN SIMULATION!  -------
res = run_simulation()

# PLOT RESULT --------
val_cols = c(
  "PLUS","MINUS",
  "MSOWL",
  "ICO",
  "COX"
)


dt.orig = setDT(data.frame(cres))

dt = copy(dt.orig)
dt[, (names(dt)) := lapply(.SD, as.numeric)]
dt[COMP_METHOD == 1, div.method := "DIV" ]
dt[COMP_METHOD == 2, div.method := "PRE" ]


# NA COUNT
ico.na.rate   = paste0(round(mean(is.na(dt$ICO)), 2)*100, "%")
ico.na.rate
msowl.na.rate = paste0(round(mean(is.na(dt$MSOWL)), 2)*100, "%")
msowl.na.rate




dt.summary = setnames(dt[, 
                         sapply(.SD, 
                                function(x) list(
                                  "m"=round(mean(x, na.rm=T),3),
                                  "s"=round(sd(x, na.rm=T),3)
                                ) 
                         ),
                         by=.(div.method) , 
                         .SDcols = val_cols
],
c("comp.meth", sapply(val_cols, paste0, c(".mu", ".sd")))
)
dt.summary


dt.melt <- melt(
  dt, 
  measure.vars = val_cols,
  id.vars= c("n","scen","cens_param","pct.cens","b", "div.method")
)
dt.melt

dt.melt$model_method = interaction(dt.melt$div.method, dt.melt$variable)

dt.melt$n[1]{}

ggplot(
  dt.melt, 
  aes(x=model_method, y=value, fill=variable)
) +  
  geom_boxplot(
    # position = position_dodge
  ) +
  theme_minimal() +
  theme(
    legend.position="bottom"
  ) + ggtitle(
    paste0(
      "B=100;",
      "; n.train=",dt.melt$n[1],
      "; ico.na.rate=", ico.na.rate,
      "; msowl.na.rate=",msowl.na.rate,
      "; Censoring Rate: ", pretty_percent(dt.melt$pct.cens)
    )
  )
# SAVE TO FOLDER FOR LATEX REPORT
fname = paste0("./documents/thesis_paper/figures/div.conq.b=100.pdf")
print(paste("saving plot to file: ", fname))
ggsave(fname, width=14, height = 10)
