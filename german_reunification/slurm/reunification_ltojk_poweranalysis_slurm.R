### Clear Environment ###
rm(list = ls(all = TRUE))

### In this example: N = 17, T_0 = 31, T = 44


### Load in the R Libraries ###
library(foreach)
library(doParallel)

source("reunification_helper.R")
countries = read.csv('german_reunification.csv')

# The following was applied to the above data, updating the data.
### Fill missing data: will backfill for same country, for same variable
#library(dplyr)
#library(tidyr)

#countries  <- countries  %>% group_by(code) %>%
#  fill(c(infrate, trade, schooling, invest60, invest70,invest80, industry), .direction = "up") %>%
#  fill(c(infrate, trade, schooling, invest60, invest70,invest80, industry), .direction = "down")
#countries  <- data.frame(countries)

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
size.resampled.dataset = 14
alpha = 0.05
tau = -50800
mc.samples = 2

size.resampled.dataset = as.integer(arguments[1])
alpha = as.numeric(arguments[2])
tau = as.numeric(arguments[3])
mc.samples = as.integer(arguments[4])

### calculate c: shift for LTOJK ###
type.I.error <- function(alpha,c){
  n = size.resampled.dataset
  Q = (1-alpha)*(n-1)*(n-2) + c*(n-1)
  f = 0.5*(3 - 3/n - sqrt(9*(1 - 1/n)^2 - 12*(1 - 2/n + 2/(3*n^2)- Q/n^2 ) ) )
  floor(n*f)/n
}

cs = seq(-1,0,0.01)
errors = rep(0,length(cs))
for (i in 1:length(cs)){
  errors[i] = type.I.error(alpha, cs[i])
}


lto.shift = cs[which(errors == type.I.error(alpha,0) )[1]]+0.01

######################################################################################
######################################################################################
######################################################################################


### Setup DoParallel ###
cores = detectCores()
cl = makeCluster(cores[1]-2) #not to overload your computer
registerDoParallel(cl)

country.code = unique(countries$code)[-7] # exclude West Germany

### Run MC simulation in parallel ###
simulation_results <- foreach(mc.run = 1:mc.samples,
                              .packages = c('optimx','rgenoud','LowRankQP')) %dopar% {
                                set.seed(mc.run)
                                
                                coverages.per.dataset = matrix(0,nrow = size.resampled.dataset, ncol = time.periods)
                                ci.lengths.per.dataset = matrix(0,nrow = size.resampled.dataset, ncol = time.periods)
                                p.values.per.dataset = matrix(0,nrow = size.resampled.dataset, ncol = time.periods)
                                
                                power.per.dataset = matrix(0,nrow = size.resampled.dataset, ncol = time.periods)
                                RMSPE.pvalue.dataset = rep(0,size.resampled.dataset)
                                
                                p.values.minusc.per.dataset = matrix(0,nrow = size.resampled.dataset, ncol = time.periods)
                                power.minusc.per.dataset = matrix(0,nrow = size.resampled.dataset, ncol = time.periods)
                                RMSPE.minusc.pvalue.dataset = rep(0,size.resampled.dataset)
                                
                                new.index = sample(country.code, size.resampled.dataset, replace = FALSE) 
                                
                                for (tu.index in 1:length(new.index)) {
                                  
                                  ### resample new dataset ###
                                  
                                  new.countries.data = countries[countries$code %in% new.index,]
                                  
                                  ### this will be the new treated unit. ###
                                  treatment.unit = new.index[tu.index]
                                  
                                  ### update the dataset with fake treatment ###
                                  true.values = new.countries.data[new.countries.data$code == treatment.unit,
                                                                 c('gdp','year')]
                                  true.treatment.values = rep(tau,time.periods)
                                  mask = (new.countries.data$code == treatment.unit) & (new.countries.data$year >= 1990)
                                  new.countries.data[mask, 'gdp'] = 
                                    new.countries.data[mask, 'gdp'] + tau # just for the evaluation of power
                                  
                                  new.index.minus.treated = new.index[! new.index %in% treatment.unit]
                                  
                                  ### R_{i,j}^LTO ###
                                  nC2 = length(new.index.minus.treated)*(length(new.index.minus.treated) - 1)/2
                                  LTO.residuals = matrix(0, nrow = nC2, ncol = time.periods)
                                  LTO.RMSPEs = rep(0,nC2)
                                  
                                  ### hat(mu)_{-i,j}(treated) ###
                                  LTO.regression.on.newpoint = matrix(0, nrow = nC2, ncol = time.periods)
                                  LTO.residual.on.newpoint = matrix(0, nrow = nC2, ncol = time.periods)
                                  LTO.RMSPEs.on.newpoint = rep(0,nC2)
                                  
                                  ### Calculate Leave Two Out points
                                  pair.index = 1
                                  for (i in 1:(length(new.index.minus.treated)-1) ) {
                                    for (j in (i+1):length(new.index.minus.treated)) {
                                      
                                      jackknife.lto.units = c(new.index.minus.treated[i],new.index.minus.treated[j])
                                      to.remove = jackknife.lto.units
                                      controls = new.index.minus.treated[! new.index.minus.treated %in% to.remove]
                                      
                                      ### get SC on LTO unit ###
                                      dataprep.out = reunification.sc.dataprep(new.countries.data, jackknife.lto.units[1], controls)
                                      synth.out <- synth(data.prep.obj = dataprep.out, method = "BFGS",quadopt = "LowRankQP")
                                      res1 <- dataprep.out$Y1plot - 
                                        (dataprep.out$Y0plot %*% synth.out$solution.w)
                                      RMSPE1 = mean(res1[(T0+1):time.periods]^2)/mean(res1[1:T0]^2)
                                      
                                      ### get SC on other LTO unit ###
                                      dataprep.out = reunification.sc.dataprep(new.countries.data, jackknife.lto.units[2], controls)
                                      synth.out <- synth(data.prep.obj = dataprep.out, method = "BFGS",quadopt = "LowRankQP")
                                      res2 <- dataprep.out$Y1plot - 
                                        (dataprep.out$Y0plot %*% synth.out$solution.w)
                                      RMSPE2 =mean(res2[(T0+1):time.periods]^2)/ mean(res2[1:T0]^2)
                                      
                                      LTO.residuals[pair.index,] = pmax(abs(res1),abs(res2))
                                      LTO.RMSPEs[pair.index] = max(RMSPE1, RMSPE2)
                                      
                                      ### run SC method on the actual treatment unit ###
                                      dataprep.out = reunification.sc.dataprep(new.countries.data, treatment.unit, controls)
                                      
                                      synth.out = synth(data.prep.obj = dataprep.out, method = "BFGS",quadopt = "LowRankQP")
                                      res = dataprep.out$Y1plot - 
                                        (dataprep.out$Y0plot %*% synth.out$solution.w)
                                      LTO.residual.on.newpoint[pair.index,] = abs(res)
                                      
                                      LTO.RMSPEs.on.newpoint[pair.index] = mean(res[(T0+1):time.periods]^2)/mean(res[1:T0]^2)
                                      LTO.regression.on.newpoint[pair.index,] = (dataprep.out$Y0plot %*% synth.out$solution.w) 
                                      
                                      pair.index = pair.index +1
                                    }
                                  }
                                  
                                  ### analyze data to create confidence intervals, p values, coverage, power ###
                                  ci.upper = rep(0, time.periods)
                                  ci.lower = rep(0, time.periods)
                                  ci.lengths.over.time = rep(0, time.periods)
                                  p.values.over.time = rep(0, time.periods)
                                  p.values.minusc.over.time = rep(0, time.periods)
                                  
                                  for (t in 1:time.periods) {
                                    ### TODO should we take the rounded off quantile? ###
                                    ci.upper[t] = quantile( LTO.regression.on.newpoint[,t] + LTO.residuals[,t], 1-alpha, type = 1)[[1]]
                                    ci.lower[t] = quantile( LTO.regression.on.newpoint[,t] - LTO.residuals[,t], alpha, type = 1)[[1]]
                                    ci.lengths.over.time[t] = ci.upper[t] - ci.lower[t]
                                    
                                    p.values.over.time[t] = sum(LTO.residual.on.newpoint[,t] <= LTO.residuals[,t])/nC2
                                    p.values.minusc.over.time[t] = p.values.over.time[t] + lto.shift/(length(new.index.minus.treated) - 1)
                                  }
                                  
                                  coverages.over.time = as.integer((true.values$gdp <= ci.upper) & (true.values$gdp >= ci.lower))
                                  coverages.per.dataset[tu.index, ] = coverages.over.time
                                  ci.lengths.per.dataset[tu.index, ] = ci.lengths.over.time
                                  p.values.per.dataset[tu.index, ] = p.values.over.time
                                  
                                  RMSPE.pvalue.dataset[tu.index] = sum( LTO.RMSPEs >= LTO.RMSPEs.on.newpoint)/nC2 
                                  power.per.dataset[tu.index,] = as.integer(p.values.per.dataset[tu.index,] <= alpha)
                                  
                                  p.values.minusc.per.dataset[tu.index, ] = p.values.minusc.over.time
                                  RMSPE.minusc.pvalue.dataset[tu.index] = RMSPE.pvalue.dataset[tu.index]+ lto.shift/(length(new.index.minus.treated) - 1)
                                  power.minusc.per.dataset[tu.index,] = as.integer(p.values.minusc.per.dataset[tu.index,] <= alpha)
                                }
                                
                                ### summarize data over the resamples of the original smoking data. ###
                                power.over.time = colMeans(power.per.dataset)
                                RMSPE.power = mean(RMSPE.pvalue.dataset < alpha)
                                coverages = colMeans(coverages.per.dataset)
                                ci.lengths = colMeans(ci.lengths.per.dataset)
                                
                                power.minusc.over.time = colMeans(power.minusc.per.dataset)
                                RMSPE.minusc.power = mean(RMSPE.minusc.pvalue.dataset < alpha)
                                
                                res = list(coverages, ci.lengths, power.over.time, RMSPE.power,
                                           power.minusc.over.time,RMSPE.minusc.power, p.values.per.dataset)
                                names(res) = c("coverages","ci.lengths","power.over.time","RMSPE.power",
                                               "power.minusc.over.time","RMSPE.minusc.power","p.values.per.dataset")
                                res
                              }

stopCluster(cl)



######################################################################################
######################################################################################
######################################################################################
### save the results ###

sim_results_file = sprintf("~/Results/SyntheticControls/reunification_power_analysis_LTO_jackknife_tau%.1f_alpha%.2f_N%d.RData",
                           tau, alpha, size.resampled.dataset)
save(simulation_results, file = sim_results_file)
