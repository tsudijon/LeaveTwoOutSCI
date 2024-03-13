### Clear Environment ###
rm(list = ls(all = TRUE))

### Load in the R Libraries ###
library(foreach)
library(doParallel)

source("smoking_helper.R")
smoking = read.csv('californiaprop99.csv')[,-1]
smoking$state = as.character(smoking$state)

### Take in arguments ###
arguments <- commandArgs(trailingOnly=TRUE)
print(arguments)

### Ensure the proper arguments are inputted
if (length(arguments) != 4) {
  stop(sprintf("At least one argument must be supplied, there are %s args",length(arguments)), call.=FALSE)
}

######################################################################################
######################################################################################
######################################################################################

time.periods = 2000 - 1970 + 1
T0 = 1988 - 1970 + 1 #onset of treatment

size.resampled.dataset = 8 # these params are to help debug
alpha = 0.15
tau = -150
mc.samples = 2

size.resampled.dataset = as.integer(arguments[1])
alpha = as.numeric(arguments[2])
tau = as.numeric(arguments[3])
mc.samples = as.integer(arguments[4])

######################################################################################
######################################################################################
######################################################################################

### Setup DoParallel ###
cores = detectCores()
cl = makeCluster(cores[1]-3) #not to overload your computer
registerDoParallel(cl)

states = c(1:2,4:38)

### Run MC simulation in parallel ###
simulation_results <- foreach(mc.run = 1:mc.samples,
                              .packages = c('optimx','rgenoud','LowRankQP')) %dopar% {
                                set.seed(mc.run)                              
                                
                                p.values.per.dataset = matrix(0,nrow = size.resampled.dataset, ncol = time.periods)
                                power.per.dataset = matrix(0,nrow = size.resampled.dataset, ncol = time.periods)
                                RMSPE.pvalue.dataset = rep(0,size.resampled.dataset)
                              
                                
                                ### resample new dataset ###
                                new.index = sample(states, size.resampled.dataset, replace = FALSE) 
                                
                                for (tu.index in 1:length(new.index)) {
                                  ### create new copy of data
                                  new.smoking.data <- smoking
                                  new.smoking.data <- new.smoking.data[new.smoking.data$stateno %in% new.index,]
                                  
                                  ### this will be the new treated unit. ###
                                  treatment.unit = new.index[tu.index]   
                                  
                                  ### update the dataset with fake treatment ###
                                  true.values = new.smoking.data[new.smoking.data$stateno == treatment.unit,
                                                                 c('cigsale','year')]
                                  true.treatment.values = rep(tau,time.periods)
                                  mask = (new.smoking.data$stateno == treatment.unit) & (new.smoking.data$year >= 1989)
                                  new.smoking.data[mask, 'cigsale'] = 
                                    new.smoking.data[mask, 'cigsale'] + tau # just for the evaluation of power
                                  
                                  new.index.minus.treated = new.index[! new.index %in% treatment.unit]
                                  
                                  ### R_i^LOO. We count only the non treated units for the LOO procedure. ###
                                  LOO.RMSPEs = rep(0,length(new.index.minus.treated))
                                  Treated.RMSPEs = rep(0,length(new.index.minus.treated))
                                  
                                  ### Calculate LOO residuals
                                  for (i in 1:length(new.index.minus.treated)) {
                                    loo.unit = new.index.minus.treated[i]
                                    
                                    # exclude original treated unit in control pool.
                                    controls = new.index[! new.index %in% c(loo.unit, treatment.unit)] 
                                    
                                    ### get SC on LOO unit ###
                                    dataprep.out = smoking.sc.dataprep(new.smoking.data, loo.unit, controls)
                                    synth.out = synth(data.prep.obj = dataprep.out, method = "BFGS",quadopt = "LowRankQP")
                                    LOO.residual = dataprep.out$Y1plot - 
                                      (dataprep.out$Y0plot %*% synth.out$solution.w) # residuals of the SC on the LOO unit
                                    
                                    LOO.RMSPEs[i] = mean(LOO.residual[(T0+1):time.periods]^2)/mean(LOO.residual[1:T0]^2)
                                    
                                    ### get SC on the original treated unit
                                    dataprep.out = smoking.sc.dataprep(new.smoking.data, treatment.unit, controls)
                                    synth.out = synth(data.prep.obj = dataprep.out, method = "BFGS",quadopt = "LowRankQP")
                                    Treated.residual = dataprep.out$Y1plot - 
                                      (dataprep.out$Y0plot %*% synth.out$solution.w)
                                    
                                    Treated.RMSPEs[i] = mean(Treated.residual[(T0+1):time.periods]^2)/mean(Treated.residual[1:T0]^2)
                                  }
                                  
                                  RMSPE.pvalue.dataset[tu.index] = sum(LOO.RMSPEs > Treated.RMSPEs)/length(new.index.minus.treated)
                                  
                                }
                                
                                ### summarize data over the resamples of the original Basque data. ###
                                RMSPE.power = mean(RMSPE.pvalue.dataset < alpha)
                                
                                res = list(RMSPE.power, p.values.per.dataset)
                                names(res) = c("RMSPE.power","p.values.per.dataset")
                                res
                                
                              }
stopCluster(cl)

######################################################################################
######################################################################################
######################################################################################
### save the results ###

sim_results_file = sprintf("Results/SyntheticControls/smoking_power_analysis_loo_test_tau%.1f_alpha%.2f_N%d.RData",
                           tau, alpha, size.resampled.dataset)
save(simulation_results, file = sim_results_file)
