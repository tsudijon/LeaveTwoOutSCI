###############
#### Setup ####
###############
library(foreach)
library(doParallel)
library(Synth)


setwd("~/projects/Research/statistics_research/synthetic controls/software/german_reunification")
source("reunification_helper.R")
german.data <- read.csv("german_reunification.csv")

# West Germany is country code 7
codes <- unique(german.data$code) 
countries <- codes[!codes %in% 7]
treatment.unit <- 7

### Set the parameters for the analysis ###
N = 17
time.periods <- 2003 - 1960 + 1
T0 <- 1990 - 1960 + 1 #onset of treatment

### Fill missing data: will backfill for same country, for same variable
library(dplyr)
library(tidyr)


german.data <- german.data %>% group_by(code) %>%
  fill(c(infrate, trade, schooling, invest60, invest70,invest80, industry), .direction = "up") %>%
  fill(c(infrate, trade, schooling, invest60, invest70,invest80, industry), .direction = "down")
german.data <- data.frame(german.data)

######################################
#### Test SC works  ####
######################################

dataprep.out <- reunification.sc.dataprep(german.data, 7, countries)
synth.out <- synth(data.prep.obj = dataprep.out, method = "BFGS",quadopt = "LowRankQP")
res1 <- dataprep.out$Y1plot - 
  (dataprep.out$Y0plot %*% synth.out$solution.w)

######################################
#### Calculate Synthetic Controls ####
######################################

### Compute all synthetic controls for West Germany, leaving out two at a time.
### this is 16C2 times 3 synthetic control calculations to run.


### R_{i,j}^LTO ###
nC2 <- 16*15/2  # technically this is N-1 Choose 2
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
                              .packages = c('optimx','rgenoud','LowRankQP')) %dopar% {
                                
                                i <- LTO.indices[pair.index,1]
                                j <- LTO.indices[pair.index,2]
                                
                                jackknife.lto.units <- c(countries[i],countries[j])
                                controls <- countries[! countries %in% jackknife.lto.units]
                                
                                ### get SC on LTO unit ###
                                dataprep.out <- reunification.sc.dataprep(german.data, jackknife.lto.units[1], controls)
                                synth.out <- synth(data.prep.obj = dataprep.out, method = "BFGS",quadopt = "LowRankQP")
                                res1 <- dataprep.out$Y1plot - 
                                  (dataprep.out$Y0plot %*% synth.out$solution.w)
                                RMSPE1 <- mean(res1[(T0+1):time.periods]^2)/mean(res1[1:T0]^2)
                                
                                ### get SC on other LTO unit ###
                                dataprep.out <- reunification.sc.dataprep(german.data, jackknife.lto.units[2], controls)
                                synth.out <- synth(data.prep.obj = dataprep.out, method = "BFGS",quadopt = "LowRankQP")
                                res2 <- dataprep.out$Y1plot - 
                                  (dataprep.out$Y0plot %*% synth.out$solution.w)
                                RMSPE2 <- mean(res2[(T0+1):time.periods]^2)/ mean(res2[1:T0]^2)
                                
                                LTO.residuals <- pmax(abs(res1),abs(res2))
                                LTO.RMSPEs <- max(RMSPE1, RMSPE2)
                                
                                ### run SC method on the actual treatment unit ###
                                dataprep.out <- reunification.sc.dataprep(german.data, treatment.unit, controls)
                                
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

#saveRDS(simulation.results, "GermanReunificationLTOResults.rds")

#####################################
### Set up placebo calculation. #####
#####################################

### run SC method on the actual treatment unit ###
dataprep.out = reunification.sc.dataprep(german.data, treatment.unit, countries)
synth.out = synth(data.prep.obj = dataprep.out, method = "BFGS",quadopt = "LowRankQP")

observed.value.over.time = dataprep.out$Y1plot - 
  (dataprep.out$Y0plot %*% synth.out$solution.w)

observed.rmspe = mean(observed.value.over.time[(T0+1):time.periods]^2)/mean(observed.value.over.time[1:T0]^2)

### Setup DoParallel ###
cores <- detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer
registerDoParallel(cl)


placebo.simulation.results <- foreach(i = 1:length(countries),
                        .packages = c('optimx','rgenoud','LowRankQP')) %dopar% {
                          
                          
                          ### Calculate Leave One Out points
                          placebo.unit = countries[i]
                          controls = codes[! codes %in% placebo.unit] # included treated unit in controls when calculating other RMSPE
                          
                          ### get SC on placebo unit ###
                          dataprep.out = reunification.sc.dataprep(german.data, placebo.unit, controls)
                          synth.out = synth(data.prep.obj = dataprep.out, method = "BFGS",quadopt = "LowRankQP")
                          placebo.residuals = dataprep.out$Y1plot - 
                            (dataprep.out$Y0plot %*% synth.out$solution.w)
                          
                          placebo.RMSPE = mean(placebo.residuals[(T0+1):time.periods]^2)/mean(placebo.residuals[1:T0]^2)
                          placebo.RMSPE
                        }


stopCluster(cl)

#saveRDS(placebo.simulation.results, "GermanReunificationPlaceboResults.rds")

######################################
#### Compute indicator matrix for SA #
######################################

simulation.results <- readRDS("GermanReunificationLTOResults.rds")


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
LTO.pvalue # pvalue is 0.04166, which is less than 1/17.

# Let us declare the p-value significant at the alpha = 0.05 level.


LTO.valid.pvalue = sum(win.matrix)/((N-1)*(N-1)) + 1/(N-1)
LTO.valid.pvalue # pvalue is 0.1015.


#####################################
#### Calculate placebo p-value ######
#####################################

placebo.simulation.results <- readRDS("GermanReunificationPlaceboResults.rds")
placebo.pvalue = (sum( placebo.simulation.results > observed.rmspe) + 1)/N
# placebo pvalue is exactly 0.0588 = 1/N.


###############################################
#### Optimize weighted p-values; product ######
###############################################
library(quadprog)

alpha = 0.05

## For a range of Gamma, optimize product p-value
Gamma = 1.5
N = 17

### create G
step1 = matrix(0, nrow = N, ncol = N-1)
step1[-treatment.unit,] <- win.matrix
G = matrix(0,nrow = N, ncol = N)
G[,-treatment.unit] <- step1


e.II = matrix(0,N,N); e.II[treatment.unit,treatment.unit] = 1
A = G + alpha*diag(N) - 2*alpha*e.II

e.I = rep(0,N); e.I[treatment.unit] = 1
b = 2*alpha*e.I

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


# feed into quadprog
#solve.QP(Dmat = 2*A, dvec = b, Amat = t(Amat), bvec = bvec)

# well, quadprog doesn't support nonconvex quadratic programs.

library(gurobi)

### Do the sensitivity analysis
## For a range of Gamma, optimize product p-value
###

gammas = seq(1,1.5,0.01)

p.values.greaterthanalpha = rep(0,length(gammas))

for (i in 1:length(gammas) ) {
  Gamma = gammas[i]
  bvec = rep(c(1/(Gamma*N),-Gamma/N),N) 
  bvec = c(bvec,1,-1)
  
  ### plug into lpSolve ###
  model <- list()
  
  model$obj        <- b # linear objective
  model$Q          <- A # quadratic objective
  model$modelsense <- "max"
  model$rhs        <- bvec # constraint values
  model$sense      <- rep(">") # constraint inequalities
  model$vtype      <- "C" # continuous variables.
  model$A          <- Amat # constraint matrix
  
  # solve the optimization problem using Gurobi
  result <- gurobi(model, list(NonConvex = 2))
  
  # check if greater than alpha
  print(result$objval)
  p.values.greaterthanalpha[i] = (result$objval > alpha)
}

p.values.greaterthanalpha

idx = tail(which(p.values.greaterthanalpha == 0),1)
gammas[idx]

# 1.08 seems to be the sensitivity


###########################################
#### Optimize weighted p-values; sum ######
###########################################

library(lpSolve)

N = 17

### create G
step1 = matrix(0, nrow = N, ncol = N-1)
step1[-treatment.unit,] <- win.matrix
G = matrix(0,nrow = N, ncol = N)
G[,-treatment.unit] <- step1


E = matrix(1, nrow = N, ncol = N)
E[treatment.unit,] = 0
diag(E) = 0

e.I = rep(0,N); e.I[treatment.unit] = 1
b = 2*alpha*(N-2)*e.I

### objective vector ###
w = diag((G + t(G))%*% E) + b


alpha = 0.05

### constraints. Encodes box constraints. Last two lines encode sum pi = 1. ###
Amat = matrix(0,nrow = 2*(N+1),N)

for (i in 1:N){
  Amat[2*i - 1,i] = 1
  Amat[2*i,i] = -1
}
Amat[2*N+1,] = rep(1,N)
Amat[2*N+2,] = rep(-1,N)

Gamma = 1.15

bvec = rep(c(1/(Gamma*N),-Gamma/N),N) 
bvec = c(bvec,1,-1)

### plug into lpSolve ###
f.obj = w
f.con = Amat
f.dir = rep(">=",dim(Amat)[1])
f.rhs = bvec

lp = lp("max", f.obj, f.con, f.dir, f.rhs)

lp$objval > 2*alpha*(N-2)

### Do the sensitivity analysis
## For a range of Gamma, optimize product p-value
###

gammas = seq(1,1.5,0.01)

p.values.greaterthanalpha = rep(0,length(gammas))

for (i in 1:length(gammas) ) {
  Gamma = gammas[i]
  bvec = rep(c(1/(Gamma*N),-Gamma/N),N) 
  bvec = c(bvec,1,-1)
  
  ### plug into lpSolve ###
  f.obj = w
  f.con = Amat
  f.dir = rep(">=",dim(Amat)[1])
  f.rhs = bvec
  
  lp = lp("max", f.obj, f.con, f.dir, f.rhs)
  
  p.values.greaterthanalpha[i] = (lp$objval > 2*alpha*(N-2))
}

idx = tail(which(p.values.greaterthanalpha == 0),1)
gammas[idx]

# 1.19 seems to be the sensitivity

###########################################
#### Optimize weighted p-values; sum Alpha GEQ 1/N ######
###########################################

alpha = 0.1
N = 17

### create G
step1 = matrix(0, nrow = N, ncol = N-1)
step1[-treatment.unit,] <- win.matrix
G = matrix(0,nrow = N, ncol = N)
G[,-treatment.unit] <- step1


E = matrix(1, nrow = N, ncol = N)
E[treatment.unit,] = 0
diag(E) = 0

e.I = rep(0,N); e.I[treatment.unit] = 1
b = (2*alpha*(N-1)-2)*e.I

### objective vector ###
w = diag((G + t(G))%*% E) + b


### constraints. Encodes box constraints. Last two lines encode sum pi = 1. ###
Amat = matrix(0,nrow = 2*(N+1),N)

for (i in 1:N){
  Amat[2*i - 1,i] = 1
  Amat[2*i,i] = -1
}
Amat[2*N+1,] = rep(1,N)
Amat[2*N+2,] = rep(-1,N)



gammas = seq(1,5,0.05)

p.values.lessthanalpha = rep(0,length(gammas))

for (i in 1:length(gammas) ) {
  Gamma = gammas[i]
  bvec = rep(c(1/(Gamma*N),-Gamma/N),N) 
  bvec = c(bvec,1,-1)
  
  ### plug into lpSolve ###
  f.obj = w
  f.con = Amat
  f.dir = rep(">=",dim(Amat)[1])
  f.rhs = bvec
  
  lp = lp("max", f.obj, f.con, f.dir, f.rhs)
  p.values.lessthanalpha[i] = (lp$objval > 2*alpha*(N-1) - 2)
}

p.values.lessthanalpha

idx = tail(which(p.values.greaterthanalpha == 0),1)
gammas[idx]



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
  model$modelsense <- "max"
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
  ylim(0,0.2) +
  xlab("Gamma") + 
  ylab("Maximum P-Value") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") + 
  theme_bw() + 
  coord_fixed(4)


