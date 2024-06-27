#============================================
#Testing Thorson's 2024 AIC with EFD script
#implementing cAIC from Thorson 2024 paper
#https://pubmed.ncbi.nlm.nih.gov/38859712/
#Catarina wor
#August 2022
#============================================
library(TMB)
source("calculate_EDF.R")


TMB::compile('Ricker_thorson.cpp')
dyn.load(dynlib('Ricker_thorson'))

#load in simulated data
thdat<-list()
thdat$obs_S = c(12505.126, 20111.116, 172686.610, 71920.143, 12651.322, 21015.153,
  23944.818, 214044.743, 5074.376,  40162.055,  74313.534,  31235.376,  12853.933,
  52286.585,  29966.241, 141410.135, 48743.597,  45806.243,  38195.634, 119869.269,
 130761.487, 138193.814, 444722.064, 191114.100,  33495.952, 124545.345,  65722.640,  
 58798.785,  70501.600,  79221.207,  42219.200,  27949.292, 57144.641,  10230.104,
   42367.015,  62051.647,   8937.645,  21204.841,  51388.608,  65003.894)
)
thdat$obs_logRS = c(1.2221039, 1.3015935, -0.2387457,  1.2111039,  0.5124565,
  1.4412325,  1.3564489, -0.6392406,  1.4492694,  0.8206990,  0.1459697,  1.8683927,
  1.6853088,  1.2074952,  0.5226575,  0.1549174,  1.5302777 , 1.7502084,  2.8313796,
  -0.3333775,  1.0137373,  0.2105126, -1.2668415, -0.1767572,  0.4675784,  0.4179322,
  0.6673385,  0.4570888, 0.6753502, -0.2781308,  1.0776220,  0.9316289,  0.6517002,
  0.1547829, -0.3928560,  1.4001415,  1.1375043,  1.1507442,  1.6665192,  1.3641949)
thdat$options_z<-c(1,0)


thparams <- list( "alpha" = 0,
                 "beta" = 0,
                 "ln_sigA" = 0,
                 "ln_sigB" = 0,
                 "ln_sigma" = 0,
                 "epsA_t" = rep(0,length(thdat$obs_S)),
                 "epsB_t" = rep(0,length(thdat$obs_S)) )
Map <- list()
Map$ln_sigB = factor(NA)
Map$epsB_t = factor( rep(NA,length(thparams$epsB_t)) )




objth <- TMB::MakeADFun(data = thdat, 
      parameters = thparams, 
      map = Map,
      random = c( "epsA_t"), 
      DLL = "Ricker_thorson")

 
#objth$env$beSilent()
nonvariance_fixed_effects = c("alpha","beta")

  # Optimize first time
optth <- nlminb( start=objth$par, obj=objth$fn, grad=objth$gr )
SD = sdreport( objth)
fit10 = list( "Opt"=optth,
              "SD" = SD,
              "Obj"=objth,
              "nonvariance_fixed_effects"=nonvariance_fixed_effects )

EDF10 = calculate_EDF(  obj=fit10$Obj,
          opt=fit10$Opt,
          nonvariance_fixed_effects=fit10$nonvariance_fixed_effects,
          prediction_name = "yhat_t",
          data_name = "obs_logRS",
          delta = 0.0001,
          show_progress = FALSE,
          refit = "full")

