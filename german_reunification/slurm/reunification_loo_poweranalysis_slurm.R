### Clear Environment ###
rm(list = ls(all = TRUE))

### In this example: N = 17, T_0 = 31, T = 44


### Load in the R Libraries ###
library(foreach)
library(doParallel)

source("reunification_helper.R")
countries = read.csv('german_reunification.csv')

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

### Set the parameters for the analysis ###
time.periods = 2003 - 1960 + 1
T0 = 1990 - 1960 + 1 #onset of treatment


# parameters below for testing purposes
# parameters below for testing purposes
size.resampled.dataset = 14
alpha = 0.05
tau = -10800
mc.samples = 2


size.resampled.dataset = as.integer(arguments[1])
alpha = as.numeric(arguments[2])
tau = as.numeric(arguments[3])
mc.samples = as.integer(arguments[4])


### calculate conditional Type I error guarantee
type.I.error <- function(alpha){
  n = size.resampled.dataset
  Q = (1-alpha)*(n-1)*(n-2)
  f = 0.5*(3 - 3/n - sqrt(9*(1 - 1/n)^2 - 12*(1 - 2/n + 2/(3*n^2)- Q/n^2 ) ) )
  floor(n*f)/n
}

LOO.alpha = (type.I.error(alpha)*size.resampled.dataset)/(2*size.resampled.dataset - 2) - 0.001 

######################################################################################
######################################################################################
######################################################################################

### Setup DoParallel ###
cores = detectCores()
cl = makeCluster(cores[1]-3) #not to overload your computer
registerDoParallel(cl)


country.code = unique(countries$code)[-7] # exclude West Germany

### Run MC simulation in parallel ###
simulation_results <- foreach(mc.run = 1:mc.samples,
                              .packages = c('optimx','rgenoud','LowRankQP')) %dopar% {
                                set.seed(mc.run)                              
                                
                                
                                p.values.per.dataset = matrix(0,nrow = size.resampled.dataset, ncol = time.periods)
                                power.per.dataset = matrix(0,nrow = size.resampled.dataset, ncol = time.periods)
                                RMSPE.pvalue.dataset = rep(0,size.resampled.dataset)
                                
                                
                                ### resample new dataset ###
                                new.index = sample(country.code, size.resampled.dataset, replace = FALSE) 
                                
                                for (tu.index in 1:length(new.index)) {
                                  ### create new copy of data
                                  new.reunification.data <- countries
                                  new.reunification.data <- new.reunification.data[new.reunification.data$code %in% new.index,]
                                  
                                  ### this will be the new treated unit. ###
                                  treatment.unit = new.index[tu.index]   
                                  
                                  ### update the dataset with fake treatment ###
                                  true.values = new.reunification.data[new.reunification.data$code == treatment.unit,
                                                                       c('gdp','year')]
                                  true.treatment.values = rep(tau,time.periods)
                                  mask = (new.reunification.data$code == treatment.unit) & (new.reunification.data$year >= 1990)
                                  new.reunification.data[mask, 'gdp'] = 
                                    new.reunification.data[mask, 'gdp'] + tau # just for the evaluation of power
                                  
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
                                    dataprep.out = reunification.sc.dataprep(new.reunification.data, loo.unit, controls)
                                    synth.out = synth(data.prep.obj = dataprep.out, method = "BFGS",quadopt = "LowRankQP")
                                    LOO.residual = dataprep.out$Y1plot - 
                                      (dataprep.out$Y0plot %*% synth.out$solution.w) # residuals of the SC on the LOO unit
                                    
                                    LOO.RMSPEs[i] = mean(LOO.residual[(T0+1):time.periods]^2)/mean(LOO.residual[1:T0]^2)
                                    
                                    ### get SC on the original treated unit
                                    dataprep.out = reunification.sc.dataprep(new.reunification.data, treatment.unit, controls)
                                    synth.out = synth(data.prep.obj = dataprep.out, method = "BFGS",quadopt = "LowRankQP")
                                    Treated.residual = dataprep.out$Y1plot - 
                                      (dataprep.out$Y0plot %*% synth.out$solution.w)
                                    
                                    Treated.RMSPEs[i] = mean(Treated.residual[(T0+1):time.periods]^2)/mean(Treated.residual[1:T0]^2)
                                  }
                                  
                                  RMSPE.pvalue.dataset[tu.index] = sum(LOO.RMSPEs > Treated.RMSPEs)/length(new.index.minus.treated)
                                  
                                }
                                
                                ### summarize data over the resamples of the original Basque data. ###
                                RMSPE.power = mean(RMSPE.pvalue.dataset < LOO.alpha)
                                
                                res = list(RMSPE.power, p.values.per.dataset)
                                names(res) = c("RMSPE.power","p.values.per.dataset")
                                res
                                
                              }

stopCluster(cl)

######################################################################################
######################################################################################
######################################################################################
### save the results ###

sim_results_file = sprintf("~/Results/SyntheticControls/reunification_power_analysis_loo_test_tau%.1f_alpha%.2f_N%d.RData",
                           tau, alpha, size.resampled.dataset)
save(simulation_results, file = sim_results_file)
