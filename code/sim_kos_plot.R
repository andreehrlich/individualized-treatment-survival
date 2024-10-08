# sim_kos_plot.R
library(ggplot2)
library(ggh4x)
library(data.table)

# PLOT RESULTS  ----------------------------------------
# Plot ITR Comparison from Simulated Data
# - compare the expected value from several itr methods, 
# - B repetitions in each (treatment fx scenario, censoring_rate, sample_size n) 
# - Boxplots 
# - latex table 
load.files.from.dir <- function(
    keyword,
    n.train.vec =  c(100,200,400),
    cens.vec = c(F,-1,-0.5,-0.2),
    scen.vec = 1:4
) {
  # LOAD FROM MULTIPLE FILES 
  ctotal = length(cens.vec)*length(scen.vec)*length(n.train.vec)
  missed = 0
  dt.kos = data.table()
  for (c_cens in cens.vec){
    for (scenario in scen.vec){
      for (n.train in n.train.vec) {
        dfile = paste0("./sim_results/", keyword, "/c", c_cens,"_s",scenario,"_n", n.train, ".csv")
        print(paste0("FILE: ", dfile))
        
        if (!file.exists(dfile)){
          print(paste("- Doesn't exist. Next."))
          missed = missed + 1
          next
        }
        
        cdt = tryCatch({
          setDT(read.csv(dfile))
        }, warning = function(w) {
          print(w)
        }, error = function(e) {
          print(e)
        }, finally = {
        })
        
        if ("error" %in% class(cdt)){
          print("- Cannot read file. Next.")
          missed = missed + 1
          next
        }
        # merge 
        if (all(dim(dt.kos) == c(0,0))){
          dt.kos = cdt
        } else{
          print(paste0(c("dt.kos",dim(dt.kos))))
          dt.kos = merge(dt.kos, cdt, all=T)
          print(dim(dt.kos))
        }
      }
    }
  }
  print(paste("missed: ", missed, "/", ctotal))
  return(dt.kos)
}

# BIG BOXPLOT ------------------------------------------------------------
plot.sim.big.boxplots <- function(dt.kos,keyword){
  
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
    # geom_line() +
    facet_wrap(~scen) +
    theme(legend.position="bottom") + 
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
  
  
  # get lambda per plot 
  lambda_df <- dt.kos[, .N, by = .(n, scen, cens.mean, lambda)]
  lambda_df
  
  ggplot(
    lambda_df,
    aes(x=log(lambda), y=N)
  ) + 
    geom_col() +
    facet_nested(
     "Sample Size" + n + "Scenario" + scen ~ "Avg Pct Censored" + cens.mean
    ) + theme() +
    ggtitle(
      paste0(
        "MSOWL; K=5; B=", max(dt.kos$b)
      ) 
    )
  fname = paste0("./documents/thesis_paper/figures/kos_big_box_lambda_", keyword, ".pdf")
  print(paste("saving plot to file: ", fname))
  ggsave(fname, width=14, height = 10)
  
  {
    # MELT data.frame to one observation per row (long format)
    id_vars = c("n", "scen", "cens_param", "pct.cens", "cens.mean","b","dur", "lambda")
    mvars = c("PLUS","MINUS", "MSOWL", "ICO", "COX") 
    dt_long <- melt(
      dt.kos, 
      measure.vars = mvars,
      id.vars= id_vars
    )
    
    # FORMAT TITLES 
    dt_long$cens.title = "Avg Pct Censored"
    dt_long[scen %in% c(1,2), surv.distr := "COX"] 
    dt_long[scen %in% c(3,4), surv.distr := "AFT"] 
    dt_long[scen %in% c(5,6,7,8), surv.distr := "EXP"] 
    dt_long[scen %in% c(9), surv.distr := "COX-NL"] 
    dt_long$surv.distr = factor(dt_long$surv.distr, levels = c("COX","AFT","EXP","COX-NL"))
    dt_long$lambda.title = "Lambda"
    # dt_long$n.label = "Sample Size"
    
    
    # BIG GRID BOXPLOT
    ggplot(
      dt_long, 
      aes(x=factor(n), y=log(value), color=variable)
    ) +  
      geom_boxplot() + 
      facet_nested(
        surv.distr + scen ~ cens.title + cens.mean,
        scales = "free"
      ) +
      theme(
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="bottom",
        legend.title=element_blank()
      ) + ggtitle(
        paste0(
          "B=", max(dt.kos$b), 
          "; ico.na.rate=", ico.na.rate,
          "; msowl.na.rate=",msowl.na.rate #,
        )
      )
    # SAVE TO FOLDER FOR LATEX REPORT
    fname = paste0("./documents/thesis_paper/figures/kos_big_box_", keyword, ".pdf")
    print(paste("saving plot to file: ", fname))
    ggsave(fname, width=14, height = 10)
  }
  
  {
    # MELT data.frame to one observation per row (long format)
    id_vars = c("n", "scen", "cens_param", "pct.cens", "cens.mean","b","dur", "lambda")
    mvars = c("MSOWL", "ICO", "COX")   # "PLUS","MINUS") 
    dt_long <- melt(
      dt.kos, 
      measure.vars = mvars,
      id.vars= id_vars
    )
    
    # FORMAT TITLES 
    dt_long$cens.title = "Avg Pct Censored"
    dt_long[scen %in% c(1,2), surv.distr := "COX"] 
    dt_long[scen %in% c(3,4), surv.distr := "AFT"] 
    dt_long[scen %in% c(5,6,7,8), surv.distr := "EXP"] 
    dt_long[scen %in% c(9), surv.distr := "COX-NL"] 
    dt_long$surv.distr = factor(dt_long$surv.distr, levels = c("COX","AFT","EXP","COX-NL"))
    dt_long$lambda.title = "Lambda"
    # dt_long$n.label = "Sample Size"
    
    
    # BIG GRID BOXPLOT
    ggplot(
      dt_long, 
      aes(x=factor(n), y=log(value), color=variable)
    ) +  
      geom_boxplot() + 
      facet_nested(
        surv.distr + scen ~ cens.title + cens.mean,
        scales = "free"
      ) +
      theme(
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="bottom",
        legend.title=element_blank()
      ) + ggtitle(
        paste0(
          "B=", max(dt.kos$b), 
          "; ico.na.rate=", ico.na.rate,
          "; msowl.na.rate=",msowl.na.rate #,
        )
      )
    # SAVE TO FOLDER FOR LATEX REPORT
    fname = paste0("./documents/thesis_paper/figures/kos_big_box_", keyword, "_no_zom.pdf")
    print(paste("saving plot to file: ", fname))
    ggsave(fname, width=14, height = 10)
  }
  
  ## LATEX TABLE ------------------------------------------------------------
  dt.kos.table <- setnames(
    dt.kos[
      ,
      sapply(.SD, function(x) list("mu"=paste0(round(mean(x,na.rm=T), 2), " (",round(sd(x,na.rm=T), 2),")"))),
      by= .(scen, cens.mean, n),
      .SDcols = mvars
    ],
    c("scen", "cens", "n",  mvars) # sapply(mvars, paste0, c(".mu")))
  )
  dt.kos.table = dt.kos.table[order(scen,cens,n)]
  
  hline_rows = c()
  for (cur.scen in 1:4) {
    res = which(dt.kos.table$scen == cur.scen)
    if (length(res) == 0) next
    hline_rows = c(hline_rows, min(res) - 1)
  }
  hline_rows
  
  # setDT(dt_mcstats)
  large <- function(x){
    paste0('{\\Large{\\bfseries ', x, '}}')
  }
  table.fname = paste0("./documents/thesis_paper/kos_table_",keyword,".tex")
  print(paste("writing xtable to file: ", table.fname))
  x = xtable(
    dt.kos.table,
    align="|r|rrl|lllll|",
    caption = c(
      "Kosorok Monte Carlo Simulation Results -- Mean (Std.Dev.) ",
      "Kosorok Simulation Results"
    )
  )
  
  # WRITE LATEXTO FILE
  print(
    x
    , include.rownames=FALSE
    , file=table.fname
    , hline.after=c(-1, 0, hline_rows, dim(dt.kos.table)[[1]])
    # , tabular.environment = "longtable"
  )
  # tbl <-ftable(dt.kos.table,row.vars=c(1,2))
  # xftbl <- xtableFtable(tbl, method = "compact")
  # print.xtableFtable(xftbl, booktabs = T)
}

# Small Boxplot ----------------------------------------------------------------
plot.sim.small.boxplots <- function(dt.kos, dt_long){
  # LATEX TABLE OF RESULTS 
  ind_cols = c(
    "scen",
    "cens.mean",
    "n"
  )
  for (cens.cur in cens.vec){
    ggplot(
      dt_long[cens_param == cens.cur, ], 
      aes(x=factor(n), y=value, color=variable)
    ) +  
      geom_boxplot() + 
      facet_wrap(~ scen, scales = "free") +
      theme(
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="bottom"
      ) + ggtitle(
        paste0(
          "B=1000;", 
          # "; n.train=",n.train,
          "; ico.na.rate=", ico.na.rate,
          "; msowl.na.rate=",msowl.na.rate,
          "; Censoring Rate: ", paste0(round(mean(dt_long$pct_cens),2)*100,"%")
        )
      )
    # SAVE TO FOLDER FOR LATEX REPORT
    fname = paste0("./documents/thesis_paper/figures/kos_sample_size_scen_value_cens_",cens.cur, "_", keyword, ".pdf")
    print(paste("saving plot to file: ", fname))
    # ggsave(fname, width=14, height = 10)
  }
}
