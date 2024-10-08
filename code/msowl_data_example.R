source("./models/msowl_fresh2.R")

# plot multi-state survival curves 
msowl_ms_plots <- function(dat, my_title){
  
  ### SURVIVAL PLOTS
  par(mfrow=c(1,2))
  # split data by live states
  dat.s1 = dat[dat$s1 == 1, ]
  dat.s1$status = factor(
    dat.s1$s2, 
    levels = 1:3, 
    labels=c("censor","response","progress/death")
  )
  
  fit <- survfit(
    Surv(time = t1, time2 = t2,event =  status) ~ A, 
    data = dat.s1,
    id=id
  )
  str(fit)
  print(fit)
  print(fit$transitions)
  
  plot(
    fit, 
    col=c(1,2,1,2), 
    lty=c(1,1,2,2),
    xlab="Days from randomisation", 
    ylab="Probability in State",
    main="State 1 (Initial Disease)"
  )
  legend(
    # x=1.8,
    # y=0.5,
    "bottomright",
    legend=c("resp:A=-1", "resp:A=+1", "death:A=-1", "death:A=+1"),
    col=c(1,2,1,2), 
    lty=c(1,1,2,2),
    lwd=0.5,
    bty="n"
  )
  
  
  # State 2 
  dat.s2 = dat[dat$s1 == 2, ]
  dat.s2$status = 1*(
    dat.s2$s1 != dat.s2$s2)
    
  # factor(
  #   dat.s2,
  #   levels = c(2,3), 
  #   labels=c("censor","progress/death")
  # )
  
  fit.s2 <- survfit(
    Surv(time = t1, time2 = t2,event =  status) ~ A, 
    data = dat.s2,
    id=id
  )
  str(fit.s2)
  print(fit.s2)
  print(fit.s2$transitions)
  
  plot(
    fit.s2, 
    col=c(1,2,1,2), 
    lty=c(1,1,2,2),
    xlab="Days since entering State", 
    ylab="Survival Probability",
    main="State 2 (Tumor Response)"
  )
  legend(
    # x=1.8,
    # y=0.5,
    "topright",
    legend=c("A=-1", "A=+1"),
    col=c(1,2),
    # lty=c(1,2,2),
    lwd=1,
    bty="n"
  )
  
  mtext(
    my_title, 
    side = 3,
    # line = -21, 
    outer = TRUE
  )
  
}


run_msowl = function(dat, cformula, cweights, tau){
  print("tune lambda parameter")
  tuned_lambda = select.lambda(
    cformula,
    id = "id",
    w = cweights,
    tau = tau,
    data = dat,
    lambdas = c(0.01, 0.1, 0.5, 1, 5, 10, 50, 100)
  )
  print(example.lambda)
  
  print("fit model")
  fit <- msowl(
    cformula,
    id = "id",
    w = cweights,
    tau = tau,
    data = dat,
    lambda = tuned_lambda$lambda
  )
  str(fit)
  print(fit)
  
  print("evaluate model")
  val <- msowl.val(
    cformula,
    id = "id",
    w = cweights,
    tau = tau,
    rule = fit,
    fixed.rules = TRUE,
    data = dat
  )
  print(val)
  
  return(list(
    lambda = tuned_lambda,
    fit = fit,
    val = val
  ))

}


# MSOWL Example Data Set 
example.data = read.csv("./models/msowl/data/example_data.csv")
# msowl_ms_plots(example.data)

# SPECTRUM Overall Survival
os.data = read.csv("./data/spectrum_ms.csv")
# msowl_ms_plots(os.data, "SPECTRUM OS")



##### ITR Analysis
w_pfs = c(1,1,0)
w_resp = c(0,1,0)
tau = 3
example.formula = ms(t1 = t1, t2 = t2, s1 = s1, s2 = s2, a = A) ~ Z1 + Z2


res.resp = run_msowl(example.data, example.formula, w_resp, tau=3)
res.pfs = run_msowl(example.data, example.formula, w_pfs, tau=3)#, tuned_lambda = 0.1)




