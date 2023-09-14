######################################################################################
######################################################################################
######################################################################################

library(Synth)

setwd("~/projects/Research/statistics_research/synthetic controls/software/basque_analysis")
source("basque_helper.R")
data("basque")

library(foreach)
library(doParallel)


######################################################################################
######################################################################################
######################################################################################

time.periods = 1997-1955+1
T0 = 1970 - 1955 + 1 #onset of treatment
mc.samples = 2
size.resampled.dataset = 8
alpha = 0.15
tau = -50

######################################################################################
######################################################################################
######################################################################################

### Setup DoParallel ###
cores = detectCores()
cl = makeCluster(cores[1]-3) #not to overload your computer
registerDoParallel(cl)

### Run MC simulation in parallel ###
simulation_results <- foreach(mc.run = 1:mc.samples,
                              .packages = c('Synth','optimx','rgenoud','kernlab')) %dopar% {
                                set.seed(mc.run)                              
                                
                                power.per.dataset = matrix(0,nrow = size.resampled.dataset, ncol = time.periods)
                                inexact.power.per.dataset = matrix(0,nrow = size.resampled.dataset, ncol = time.periods)
                                
                                p.values.per.dataset = matrix(0,nrow = size.resampled.dataset, ncol = time.periods)
                                inexact.p.values.per.dataset = matrix(0,nrow = size.resampled.dataset, ncol = time.periods)
                                
                                RMSPE.pvalue.dataset = rep(0,size.resampled.dataset)
                                inexact.RMSPE.pvalue.dataset = rep(0,size.resampled.dataset)
                                
                                ### resample new dataset ###
                                new.index = sample(c(2:16,18), size.resampled.dataset, replace = FALSE) 
                                
                                for (tu.index in 1:length(new.index)) {
                                  ### create new copy of data
                                  new.basque.data <- basque
                                  new.basque.data <- new.basque.data[new.basque.data$regionno %in% new.index,]
                                  
                                  ### this will be the new treated unit. ###
                                  treatment.unit = new.index[tu.index]   
                                  
                                  ### update the dataset with fake treatment ###
                                  true.values = new.basque.data[new.basque.data$regionno == treatment.unit,
                                                                c('gdpcap','year')]
                                  true.treatment.values = rep(tau,time.periods)
                                  mask = (new.basque.data$regionno == treatment.unit) & (new.basque.data$year >= 1970)
                                  new.basque.data[mask, 'gdpcap'] = 
                                    new.basque.data[mask, 'gdpcap'] + tau # just for the evaluation of power
                                  
                                  new.index.minus.treated = new.index[! new.index %in% treatment.unit]
                                  
                                  ### R_i^LOO. We count only the non treated units for the LOO procedure. ###
                                  LOO.residuals = matrix(0, nrow = length(new.index.minus.treated), ncol = time.periods)
                                  LOO.RMSPEs = rep(0,length(new.index.minus.treated))
                                  
                                  ### Calculate Leave One Out points
                                  for (i in 1:length(new.index.minus.treated)) {
                                    loo.unit = new.index.minus.treated[i]
                                    controls = new.index[! new.index %in% loo.unit] # included treated unit in controls when calculating other RMSPE
                                    
                                    ### get SC on LOO unit ###
                                    dataprep.out = basque.sc.dataprep(new.basque.data, loo.unit, controls)
                                    synth.out = synth(data.prep.obj = dataprep.out, method = "BFGS",quadopt = "LowRankQP")
                                    LOO.residuals[i,] = dataprep.out$Y1plot - 
                                      (dataprep.out$Y0plot %*% synth.out$solution.w)
                                    
                                    LOO.RMSPEs[i] = mean(LOO.residuals[i,(T0+1):time.periods]^2)/mean(LOO.residuals[i,1:T0]^2)
                                  }
                                  
                                  ### run SC method on the actual treatment unit ###
                                  dataprep.out = basque.sc.dataprep(new.basque.data, treatment.unit, new.index.minus.treated)
                                  synth.out = synth(data.prep.obj = dataprep.out, method = "BFGS",quadopt = "LowRankQP")
                                  
                                  observed.value.over.time = dataprep.out$Y1plot - 
                                    (dataprep.out$Y0plot %*% synth.out$solution.w)
                                  
                                  observed.rmspe = mean(observed.value.over.time[(T0+1):time.periods]^2)/mean(observed.value.over.time[1:T0]^2)
                                  
                                  ### analyze data to create p values, get power. ###
                                  for (t in 1:time.periods) {
                                    p.values.per.dataset[tu.index,t] = 
                                      (sum( abs(LOO.residuals[,t]) >= abs(observed.value.over.time[t]) ) + 1)/length(new.index)
                                    
                                    inexact.p.values.per.dataset[tu.index,t] = 
                                      (sum( abs(LOO.residuals[,t]) >= abs(observed.value.over.time[t]) ))/length(new.index)
                                  }
                                  
                                  RMSPE.pvalue.dataset[tu.index] = (sum( LOO.RMSPEs > observed.rmspe) + 1)/length(new.index)
                                  inexact.RMSPE.pvalue.dataset[tu.index] = (sum( LOO.RMSPEs > observed.rmspe))/length(new.index)
                                  
                                  power.per.dataset[tu.index ,] = as.integer(p.values.per.dataset[tu.index,] <= alpha)
                                  inexact.power.per.dataset[tu.index ,] = 
                                    as.integer(inexact.p.values.per.dataset[tu.index,] <= alpha)
                                  
                                }
                                
                                ### summarize data over the resamples of the original Basque data. ###
                                power.over.time = colMeans(power.per.dataset)
                                RMSPE.power = mean(RMSPE.pvalue.dataset <= alpha) 
                                
                                inexact.power.over.time = colMeans(inexact.power.per.dataset)
                                inexact.RMSPE.power = mean(inexact.RMSPE.pvalue.dataset <= alpha) 
                                
                                res = list(power.over.time, RMSPE.power,
                                           inexact.power.over.time, inexact.RMSPE.power)
                                names(res) = c("power.over.time","RMSPE.power",
                                               "inexact.power.over.time","inexact.RMSPE.power")
                                res
                              }

stopCluster(cl)

######################################################################################
######################################################################################
######################################################################################
### save the results ###

sim_results_file = sprintf("Results/basque_power_analysis_placebo_test_tau%.1f_alpha%.2f_N%d.RData",
                           tau, alpha, size.resampled.dataset)
save(simulation_results, file = sim_results_file)

