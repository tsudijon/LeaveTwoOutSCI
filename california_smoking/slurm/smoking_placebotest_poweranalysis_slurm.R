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
        
        power.per.dataset = matrix(0,nrow = size.resampled.dataset, ncol = time.periods)
        inexact.power.per.dataset = matrix(0,nrow = size.resampled.dataset, ncol = time.periods)
        randomized.power.per.dataset = matrix(0,nrow = size.resampled.dataset, ncol = time.periods)
        
        p.values.per.dataset = matrix(0,nrow = size.resampled.dataset, ncol = time.periods)
        inexact.p.values.per.dataset = matrix(0,nrow = size.resampled.dataset, ncol = time.periods)
        randomized.p.values.per.dataset = matrix(0,nrow = size.resampled.dataset, ncol = time.periods)
        
        RMSPE.pvalue.dataset = rep(0,size.resampled.dataset)
        inexact.RMSPE.pvalue.dataset = rep(0,size.resampled.dataset)
        randomized.RMSPE.pvalue.dataset = rep(0,size.resampled.dataset)
        
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
          LOO.residuals = matrix(0, nrow = length(new.index.minus.treated), ncol = time.periods)
          LOO.RMSPEs = rep(0,length(new.index.minus.treated))
          
          ### Calculate Leave One Out points
          for (i in 1:length(new.index.minus.treated)) {
            loo.unit = new.index.minus.treated[i]
            controls = new.index[! new.index %in% loo.unit] # included treated unit in controls when calculating other RMSPE
            
            ### get SC on LOO unit ###
            dataprep.out = smoking.sc.dataprep(new.smoking.data, loo.unit, controls)
            synth.out = synth(data.prep.obj = dataprep.out, method = "BFGS",quadopt = "LowRankQP")
            LOO.residuals[i,] = dataprep.out$Y1plot - 
              (dataprep.out$Y0plot %*% synth.out$solution.w)
            
            LOO.RMSPEs[i] = mean(LOO.residuals[i,(T0+1):time.periods]^2)/mean(LOO.residuals[i,1:T0]^2)
          }
          
          ### run SC method on the actual treatment unit ###
          dataprep.out = smoking.sc.dataprep(new.smoking.data, treatment.unit, new.index.minus.treated)
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
            
            randomized.p.values.per.dataset[tu.index,t] = 
              p.values.per.dataset[tu.index,t] - runif(1)/length(new.index)
            
          }
          
          
          ### P-value via RMSPE statistic ###
          RMSPE.pvalue.dataset[tu.index] = (sum( LOO.RMSPEs > observed.rmspe) + 1)/length(new.index)
          inexact.RMSPE.pvalue.dataset[tu.index] = (sum( LOO.RMSPEs > observed.rmspe))/length(new.index)
          randomized.RMSPE.pvalue.dataset[tu.index] = RMSPE.pvalue.dataset[tu.index] - runif(1)/length(new.index)
          
          power.per.dataset[tu.index ,] = as.integer(p.values.per.dataset[tu.index,] <= alpha)
          inexact.power.per.dataset[tu.index ,] = 
            as.integer(inexact.p.values.per.dataset[tu.index,] <= alpha)
          randomized.power.per.dataset[tu.index ,] = 
            as.integer(randomized.p.values.per.dataset[tu.index,] <= alpha)
        }
        
        ### summarize data over the resamples of the original smoking data. ###
        power.over.time = colMeans(power.per.dataset)
        RMSPE.power = mean(RMSPE.pvalue.dataset <= alpha) 
        
        inexact.power.over.time = colMeans(inexact.power.per.dataset)
        inexact.RMSPE.power = mean(inexact.RMSPE.pvalue.dataset <= alpha) 
        
        randomized.power.over.time = colMeans(randomized.power.per.dataset)
        randomized.RMSPE.power = mean(randomized.RMSPE.pvalue.dataset <= alpha) 
        
        res = list(power.over.time, RMSPE.power,
                   inexact.power.over.time, inexact.RMSPE.power,
                   randomized.power.over.time, randomized.RMSPE.power)
        names(res) = c("power.over.time","RMSPE.power",
                       "inexact.power.over.time","inexact.RMSPE.power",
                       "randomized.power.over.time","randomized.RMSPE.power")
        res
        
        }
stopCluster(cl)

######################################################################################
######################################################################################
######################################################################################
### save the results ###

sim_results_file = sprintf("~/Results/SyntheticControls/smoking_power_analysis_random_placebo_test_tau%.1f_alpha%.2f_N%d.RData",
                           tau, alpha, size.resampled.dataset)
save(simulation_results, file = sim_results_file)

