library(Synth)
source("basque_helper.R")
data("basque")

library(foreach)
library(doParallel)

cores = detectCores()
cl = makeCluster(cores[1]-3) #not to overload your computer
registerDoParallel(cl)

######################################################################################
######################################################################################
######################################################################################


### Set the parameters for the analysis ###
time.periods = 1997-1955+1

T0 = 1970 - 1955 + 1 #onset of treatment
mc.samples = 9
size.resampled.dataset = 8
alpha = 0.125
tau = 1

######################################################################################
######################################################################################
######################################################################################


### Setup DoParallel ###
cores = detectCores()
cl = makeCluster(cores[1]-3) #not to overload your computer
registerDoParallel(cl)

###
triples.in.order = matrix(0,nrow = N*(N-1)(N-2)/6, ncol = 3)
counter = 0
for (i in 1:(N-2)){
  for (j in (i+1):(N-1)){
    for (k in (j+1):N){
      counter = counter + 1
      triples.in.order[counter,] = c(i,j,k)
    }
  }
}

### Run MC simulation in parallel ###
for(mc.run in 1:mc.samples) {
        set.seed(mc.run)
        
        coverages.per.dataset = matrix(0,nrow = size.resampled.dataset, ncol = time.periods)
        ci.lengths.per.dataset = matrix(0,nrow = size.resampled.dataset, ncol = time.periods)
        p.values.per.dataset = matrix(0,nrow = size.resampled.dataset, ncol = time.periods)
        
        power.per.dataset = matrix(0,nrow = size.resampled.dataset, ncol = time.periods)
        RMSPE.pvalue.dataset = rep(0,size.resampled.dataset)
        
        ### resample new dataset ###
        new.index = sample(c(2:16,18), size.resampled.dataset, replace = FALSE) 
        new.basque.data = basque[basque$regionno %in% new.index,]
        
        
        ### go through every combination of increasing triples ###
        
        
        
        for (tu.index in 1:length(new.index)) {
          
          # this will be the new treated unit.
          treatment.unit = new.index[tu.index]
          
          # update the dataset with fake treatment
          true.values = new.basque.data[new.basque.data$regionno == treatment.unit,
                                        c('gdpcap','year')]
          true.treatment.values = rep(tau,time.periods)
          new.basque.data[new.basque.data$regionno == treatment.unit, 'gdpcap'] = 
            new.basque.data[new.basque.data$regionno == treatment.unit, 'gdpcap'] + tau # just for the evaluation of power
          
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
          for (i in 1:(length(new.index.minus.treated)-1) ) {
            for (j in (i+1):length(new.index.minus.treated)) {
              
              jackknife.lto.units = c(new.index.minus.treated[i],new.index.minus.treated[j])
              to.remove = jackknife.lto.units
              controls = new.index.minus.treated[! new.index.minus.treated %in% to.remove]
              
              ### get SC on LTO unit ###
              dataprep.out = basque.sc.dataprep(new.basque.data, jackknife.lto.units[1], controls)
              synth.out <- synth(data.prep.obj = dataprep.out, method = "BFGS")
              res1 <- dataprep.out$Y1plot - 
                (dataprep.out$Y0plot %*% synth.out$solution.w)
              RMSPE1 = mean(res1[1:T0]^2)/mean(res1[(T0+1):time.periods]^2)
              
              ### get SC on other LTO unit ###
              dataprep.out = basque.sc.dataprep(new.basque.data, jackknife.lto.units[2], controls)
              synth.out <- synth(data.prep.obj = dataprep.out, method = "BFGS")
              res2 <- dataprep.out$Y1plot - 
                (dataprep.out$Y0plot %*% synth.out$solution.w)
              RMSPE2 = mean(res2[1:T0]^2)/mean(res2[(T0+1):time.periods]^2)
              
              niC2 = (length(new.index.minus.treated)-i+1)*(length(new.index.minus.treated) - i + 2)/2
              pair.index = nC2 - niC2 + j - i
              LTO.residuals[pair.index,] = pmax(abs(res1),abs(res2))
              LTO.RMSPEs[pair.index] = max(RMSPE1, RMSPE2)
              
              ### run SC method on the actual treatment unit ###
              dataprep.out = basque.sc.dataprep(new.basque.data, treatment.unit, controls)
              
              synth.out = synth(data.prep.obj = dataprep.out, method = "BFGS")
              res = dataprep.out$Y1plot - 
                (dataprep.out$Y0plot %*% synth.out$solution.w)
              LTO.residual.on.newpoint[pair.index,] = res
              
              LTO.RMSPEs.on.newpoint[pair.index] = mean(res[1:T0]^2)/mean(res[(T0+1):time.periods]^2)
              LTO.regression.on.newpoint[pair.index,] = (dataprep.out$Y0plot %*% synth.out$solution.w) 
            }
          }
          
          ### analyze data to create confidence intervals, p values, coverage, power ###
          ci.upper = rep(0, time.periods)
          ci.lower = rep(0, time.periods)
          ci.lengths.over.time = rep(0, time.periods)
          p.values.over.time = rep(0, time.periods)
          
          for (t in 1:time.periods) {
            ### TODO should we take the rounded off quantile? ###
            ci.upper[t] = quantile( LTO.regression.on.newpoint[,t] + LTO.residuals[,t], 1-alpha, type = 1)[[1]]
            ci.lower[t] = quantile( LTO.regression.on.newpoint[,t] - LTO.residuals[,t], alpha, type = 1)[[1]]
            ci.lengths.over.time[t] = ci.upper[t] - ci.lower[t]
            
            p.values.over.time[t] = sum(abs(LTO.residual.on.newpoint[,t]) < LTO.residuals[,t])/nC2
          }
          
          coverages.per.dataset[tu.index, ] = as.integer((true.values$gdpcap < ci.upper) & (true.values$gdpcap > ci.lower))
          ci.lengths.per.dataset[tu.index, ] = ci.lengths.over.time
          p.values.per.dataset[tu.index, ] = p.values.over.time
          
          RMSPE.pvalue.dataset[tu.index] = sum( LTO.RMSPEs >= LTO.RMSPEs.on.newpoint)/nC2
          power.per.dataset[tu.index,] = as.integer(p.values.per.dataset[tu.index,] <= alpha)
          
        }
        
        ### summarize data over the resamples of the original Basque data. ###
        power.over.time = colMeans(power.per.dataset)
        RMSPE.power = mean(RMSPE.pvalue.dataset < alpha) 
        coverages = colMeans(coverages.per.dataset)
        ci.lengths = colMeans(ci.lengths.per.dataset)
        
        res = list(coverages, ci.lengths, power.over.time, RMSPE.power)
        names(res) = c("coverages","ci.lengths","power.over.time","RMSPE.power")
        res
  }

stopCluster(cl)