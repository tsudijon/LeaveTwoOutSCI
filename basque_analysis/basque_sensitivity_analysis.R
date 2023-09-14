###############
#### Setup ####
###############
library(foreach)
library(doParallel)
library(Synth)


setwd("~/projects/Research/statistics_research/synthetic controls/software/basque_analysis")
source("basque_helper.R")
data(basque)

# Basque country is country code 17
regionno <- unique(basque$regionno) 
regionno <- regionno[!regionno %in% 1]  # region 1 is all of Spain, not a region
regions <- regionno[!regionno %in% 17] 
treatment.unit <- 17

### Set the parameters for the analysis ###
N = 17
time.periods = 1997-1955+1
T0 = 1970 - 1955 + 1 #onset of treatment


######################################
#### Test SC works  ####
######################################

dataprep.out <- basque.sc.dataprep(basque, 17, regions)
synth.out <- synth(data.prep.obj = dataprep.out, method = "BFGS",quadopt = "LowRankQP")
res1 <- dataprep.out$Y1plot - 
  (dataprep.out$Y0plot %*% synth.out$solution.w)

######################################
#### Calculate Synthetic Controls ####
######################################

### Compute all synthetic controls for California, leaving out two at a time.
### this is 37C2 times 3 synthetic control calculations to run.


### R_{i,j}^LTO ###
nC2 <- (N-1)*(N-2)/2  # technically this is N-1 Choose 2
LTO.residuals <- matrix(0, nrow = nC2, ncol = time.periods)
LTO.RMSPEs <- rep(0,nC2)


#make indices

LTO.indices <- matrix(0, nrow = nC2, ncol = 2)
pair.index <- 1
for (i in 1:(N-2)) {
  for (j in (i+1):(N-1)) {
    LTO.indices[pair.index,] <- c(i,j)
    pair.index <- pair.index + 1
  }
}

### Setup DoParallel ###
cores <- detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer
registerDoParallel(cl)

simulation.results <- foreach(pair.index = 1:nC2,
                              .packages = c('optimx','rgenoud','LowRankQP','Synth')) %dopar% {
                                
                                i <- LTO.indices[pair.index,1]
                                j <- LTO.indices[pair.index,2]
                                
                                jackknife.lto.units <- c(regions[i],regions[j])
                                controls <- regions[! regions %in% jackknife.lto.units]
                                
                                ### get SC on LTO unit ###
                                dataprep.out <- basque.sc.dataprep(basque, jackknife.lto.units[1], controls)
                                synth.out <- synth(data.prep.obj = dataprep.out, method = "BFGS",quadopt = "LowRankQP")
                                res1 <- dataprep.out$Y1plot - 
                                  (dataprep.out$Y0plot %*% synth.out$solution.w)
                                RMSPE1 <- mean(res1[(T0+1):time.periods]^2)/mean(res1[1:T0]^2)
                                
                                ### get SC on other LTO unit ###
                                dataprep.out <- basque.sc.dataprep(basque, jackknife.lto.units[2], controls)
                                synth.out <- synth(data.prep.obj = dataprep.out, method = "BFGS",quadopt = "LowRankQP")
                                res2 <- dataprep.out$Y1plot - 
                                  (dataprep.out$Y0plot %*% synth.out$solution.w)
                                RMSPE2 <- mean(res2[(T0+1):time.periods]^2)/ mean(res2[1:T0]^2)
                                
                                LTO.residuals <- pmax(abs(res1),abs(res2))
                                LTO.RMSPEs <- max(RMSPE1, RMSPE2)
                                
                                ### run SC method on the actual treatment unit ###
                                dataprep.out <- basque.sc.dataprep(basque, treatment.unit, controls)
                                
                                synth.out <- synth(data.prep.obj = dataprep.out, method = "BFGS",quadopt = "LowRankQP")
                                res <- dataprep.out$Y1plot - 
                                  (dataprep.out$Y0plot %*% synth.out$solution.w)
                                
                                LTO.residual.on.newpoint <- abs(res)
                                LTO.RMSPEs.on.newpoint <- mean(res[(T0+1):time.periods]^2)/mean(res[1:T0]^2)
                                
                                res <- list(pair.index, LTO.residual.on.newpoint, LTO.RMSPEs.on.newpoint, LTO.residuals, LTO.RMSPEs)
                                names(res) <- c("pair.index", "LTO.time.res.treatment", "LTO.RMSPE.treatment",
                                                "LTO.time.max.res.controls", "LTO.RMSPE.max.controls")
                                
                                res
                              }


stopCluster(cl)

#saveRDS(simulation.results, "Prop99basqueLTOResults.rds")

#####################################
### Set up placebo calculation. #####
#####################################

### run SC method on the actual treatment unit ###
dataprep.out = basque.sc.dataprep(basque, treatment.unit, regions)
synth.out = synth(data.prep.obj = dataprep.out, method = "BFGS",quadopt = "LowRankQP")

observed.value.over.time = dataprep.out$Y1plot - 
  (dataprep.out$Y0plot %*% synth.out$solution.w)

observed.rmspe = mean(observed.value.over.time[(T0+1):time.periods]^2)/mean(observed.value.over.time[1:T0]^2)

### Setup DoParallel ###
cores <- detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer
registerDoParallel(cl)


placebo.simulation.results <- foreach(i = 1:length(regions),
                                      .packages = c('optimx','rgenoud','LowRankQP')) %dopar% {
                                        
                                        
                                        ### Calculate Leave One Out points
                                        placebo.unit = regions[i]
                                        controls = regionno[! regionno %in% placebo.unit] # included treated unit in controls when calculating other RMSPE
                                        
                                        ### get SC on placebo unit ###
                                        dataprep.out = basque.sc.dataprep(basque, placebo.unit, controls)
                                        synth.out = synth(data.prep.obj = dataprep.out, method = "BFGS",quadopt = "LowRankQP")
                                        placebo.residuals = dataprep.out$Y1plot - 
                                          (dataprep.out$Y0plot %*% synth.out$solution.w)
                                        
                                        placebo.RMSPE = mean(placebo.residuals[(T0+1):time.periods]^2)/mean(placebo.residuals[1:T0]^2)
                                        placebo.RMSPE
                                      }


stopCluster(cl)

#saveRDS(placebo.simulation.results, "Prop99basquePlaceboResults.rds")

######################################
#### Compute indicator matrix for SA #
######################################

simulation.results <- readRDS("Prop99basqueLTOResults.rds")


### Calculate the indicator of wins in three-way matches.
win.matrix <- matrix(0,nrow = 16, ncol = 16)

for (index in 1:nC2){
  i <- LTO.indices[index,1]
  j <- LTO.indices[index,2]
  
  RMSPE.treatment <- simulation.results[[index]]$LTO.RMSPE.treatment
  max.other.LTO.RMSPE <- simulation.results[[index]]$LTO.RMSPE.max.controls
  
  win.matrix[i,j] =  (RMSPE.treatment < max.other.LTO.RMSPE)
  win.matrix[j,i] = win.matrix[i,j]
}

win.matrix

#############################
#### Calculate LTO p-value ######
#############################

LTO.pvalue = sum(win.matrix)/(2*nC2)
LTO.pvalue # pvalue is 0.6666

# Let us declare the p-value insignificant at 0.05


LTO.valid.pvalue = sum(win.matrix)/((N-1)*(N-1)) + 1/(N-1)
LTO.valid.pvalue # pvalue is 0.6875


#####################################
#### Calculate placebo p-value ######
#####################################

placebo.simulation.results <- readRDS("Prop99basquePlaceboResults.rds")
placebo.pvalue = (sum( placebo.simulation.results > observed.rmspe) + 1)/N
placebo.pvalue
# placebo pvalue is exactly 0.4117

#################################################################
##### Trace out curve of p-values as a function of Gammas #######
#################################################################

# build out a binary search for this task


maximize.SA.problem.LTO <- function(Gamma, threshold) {
  N = 17
  
  ### create G
  step1 = matrix(0, nrow = N, ncol = N-1)
  step1[-treatment.unit,] <- win.matrix
  G = matrix(0,nrow = N, ncol = N)
  G[,-treatment.unit] <- step1
  
  e.II = matrix(0,N,N); e.II[treatment.unit,treatment.unit] = 1
  A = G + threshold*diag(N) - 2*threshold*e.II
  
  e.I = rep(0,N); e.I[treatment.unit] = 1
  b = 2*threshold*e.I
  
  # constraints. Encodes box constraints. Last two lines encode sum pi = 1.
  Amat = matrix(0,nrow = 2*(N+1),N)
  
  for (i in 1:N){
    Amat[2*i - 1,i] = 1
    Amat[2*i,i] = -1
  }
  Amat[2*N+1,] = rep(1,N)
  Amat[2*N+2,] = rep(-1,N)
  
  bvec = rep(c(1/(Gamma*N),-Gamma/N),N) 
  bvec = c(bvec,1,-1)
  
  bvec = rep(c(1/(Gamma*N),-Gamma/N),N) 
  bvec = c(bvec,1,-1)
  
  ### plug into Gurobi ###
  model <- list()
  
  model$obj        <- b # linear objective
  model$Q          <- A # quadratic objective
  model$modelsense <- "min"
  model$rhs        <- bvec # constraint values
  model$sense      <- rep(">") # constraint inequalities
  model$vtype      <- "C" # continuous variables.
  model$A          <- Amat # constraint matrix
  
  # solve the optimization problem using Gurobi
  result <- gurobi(model, list(NonConvex = 2))
  
  # return whether greater than threshold
  (result$objval > threshold)
}

determine.max.p.value <- function(Gamma, iter = 12){
  threshold = 0.5
  for (i in 1:iter) {
    greater.than.threshold = maximize.SA.problem.LTO(Gamma, threshold)
    
    if (greater.than.threshold) {
      threshold = threshold + 0.5^(i+1)
    } else {
      threshold = threshold - 0.5^i + 0.5^(i+1)
    }
  }  
  threshold
}

### Trace curve of p-values

gammas = seq(1,2,0.01)
SA.p.values = rep(0,length(gammas))

for (i in 1:length(gammas) ) {
  Gamma = gammas[i]
  
  SA.p.values[i] = determine.max.p.value(Gamma)
}
SA.p.values

library(ggplot2)

df <- data.frame(gammas = gammas, pvals = SA.p.values)
ggplot(data=df, aes(x = gammas, y = pvals)) + 
  geom_line(linewidth = 1.)+
  #geom_point() + 
  #ylim(0.3,0.7) +
  xlab("Gamma") + 
  ylab("Minimum P-Value per Fixed Gamma") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") + 
  theme_bw() + 
  ggtitle("Basque LTO Sensitivity Analysis") + 
  theme(plot.title = element_text(hjust = 0.5))
  #coord_fixed(2) 
  
  
  
  
  
  
  
  
  
  
  