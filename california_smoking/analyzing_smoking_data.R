loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

setwd("~/Projects/Research/statistics_research/LeaveTwoOutSCI/")

size.resampled.dataset = 30
alpha = 0.02
taus = c(0,-30,-60,-90)
mc.samples = 20

######################################################################################
######################################################################################
######################################################################################
### LOO ###

powers.RMSPE.LOO =  matrix(0, length(taus),mc.samples)
power.mean.RMSPE.LOO  = rep(0,length(taus))

for (j in 1:length(taus)){
  tau = taus[j]
  data.filename <- sprintf("Results/SyntheticControls/smoking_power_analysis_loo_test_tau%.1f_alpha%.2f_N%d.RData",
                           tau, alpha, size.resampled.dataset)
  d <- loadRData(data.filename)
  
  for (i in 1:mc.samples){
    powers.RMSPE.LOO[j,i] = d[[i]]$RMSPE.power
  }
}
power.mean.RMSPE.LOO = rowMeans(powers.RMSPE.LOO)

######################################################################################
######################################################################################
######################################################################################
### LTO ###

powers.RMSPE.LTO =  matrix(0, length(taus),mc.samples)
power.mean.RMSPE.LTO  = rep(0,length(taus))
power.std.RMSPE.LTO  = rep(0,length(taus))

powers.minusc.RMSPE.LTO =  matrix(0, length(taus),mc.samples)
power.minusc.mean.RMSPE.LTO  = rep(0,length(taus))
power.minusc.std.RMSPE.LTO  = rep(0,length(taus))

powers.over.time.LTO =  array(0,dim = c(length(taus), 31,mc.samples))
power.mean.over.time.LTO = matrix(0,length(taus), 31)
power.std.over.time.LTO = matrix(0,length(taus), 31)

powers.minusc.over.time.LTO =  array(0,dim = c(length(taus), 31,mc.samples))
power.minusc.mean.over.time.LTO = matrix(0,length(taus), 31)
power.minusc.std.over.time.LTO = matrix(0,length(taus), 31)


for (j in 1:length(taus)){
  tau = taus[j]
  data.filename <- sprintf("Results/SyntheticControls/smoking_power_analysis_LTO_jackknife_tau%.1f_alpha%.2f_N%d.RData",
                           tau, alpha, size.resampled.dataset)
  d <- loadRData(data.filename)
  
  for (i in 1:mc.samples){
    powers.RMSPE.LTO[j,i] = d[[i]]$RMSPE.power
    powers.over.time.LTO[j,,i] = d[[i]]$power.over.time
    
    powers.minusc.RMSPE.LTO[j,i] = d[[i]]$RMSPE.minusc.power
    powers.minusc.over.time.LTO[j,,i] = d[[i]]$power.minusc.over.time
  }
}
power.mean.RMSPE.LTO = rowMeans(powers.RMSPE.LTO)
power.std.RMSPE.LTO = apply(powers.RMSPE.LTO, 1, sd)

power.mean.over.time.LTO = apply(powers.over.time.LTO, c(1,2), mean)
power.std.over.time.LTO = apply(powers.over.time.LTO, c(1,2), sd)

power.minusc.mean.RMSPE.LTO = rowMeans(powers.minusc.RMSPE.LTO)
power.minusc.std.RMSPE.LTO = apply(powers.minusc.RMSPE.LTO, 1, sd)

power.minusc.mean.over.time.LTO = apply(powers.minusc.over.time.LTO, c(1,2), mean)
power.minusc.std.over.time.LTO = apply(powers.minusc.over.time.LTO, c(1,2), sd)

######################################################################################
######################################################################################
######################################################################################
### Placebo ###

powers.RMSPE.Placebo =  matrix(0, length(taus),mc.samples)
power.mean.RMSPE.Placebo  = rep(0,length(taus))
power.std.RMSPE.Placebo  = rep(0,length(taus))

powers.over.time.Placebo  =  array(0,dim = c(length(taus), 31,mc.samples))
power.mean.over.time.Placebo  = matrix(0,length(taus), 31)
power.std.over.time.Placebo  = matrix(0,length(taus), 31)

inexact.powers.RMSPE.Placebo =  matrix(0, length(taus),mc.samples)
inexact.power.mean.RMSPE.Placebo  = rep(0,length(taus))
inexact.power.std.RMSPE.Placebo  = rep(0,length(taus))

inexact.powers.over.time.Placebo  =  array(0,dim = c(length(taus), 31,mc.samples))
inexact.power.mean.over.time.Placebo  = matrix(0,length(taus), 31)
inexact.power.std.over.time.Placebo  = matrix(0,length(taus), 31)

random.powers.RMSPE.Placebo =  matrix(0, length(taus),mc.samples)
random.power.mean.RMSPE.Placebo  = rep(0,length(taus))
random.power.std.RMSPE.Placebo  = rep(0,length(taus))

random.powers.over.time.Placebo  =  array(0,dim = c(length(taus), 31,mc.samples))
random.power.mean.over.time.Placebo  = matrix(0,length(taus), 31)
random.power.std.over.time.Placebo  = matrix(0,length(taus), 31)

for (j in 1:length(taus)){
  tau = taus[j]
  data.filename <- sprintf("Results/SyntheticControls/smoking_power_analysis_random_placebo_test_tau%.1f_alpha%.2f_N%d.RData",
                           tau, alpha, size.resampled.dataset)
  d <- loadRData(data.filename)
  
  for (i in 1:mc.samples){
    powers.RMSPE.Placebo[j,i] = d[[i]]$RMSPE.power
    powers.over.time.Placebo[j,,i] = d[[i]]$power.over.time
    
    inexact.powers.RMSPE.Placebo[j,i] = d[[i]]$inexact.RMSPE.power
    inexact.powers.over.time.Placebo[j,,i] = d[[i]]$inexact.power.over.time
    
    random.powers.RMSPE.Placebo[j,i] = d[[i]]$randomized.RMSPE.power
    random.powers.over.time.Placebo[j,,i] = d[[i]]$randomized.power.over.time
  }
}
power.mean.RMSPE.Placebo = rowMeans(powers.RMSPE.Placebo)
power.std.RMSPE.Placebo = apply(powers.RMSPE.Placebo, 1, sd)

power.mean.over.time.Placebo = apply(powers.over.time.Placebo, c(1,2), mean)
power.std.over.time.Placebo = apply(powers.over.time.Placebo, c(1,2), sd)

inexact.power.mean.RMSPE.Placebo = rowMeans(inexact.powers.RMSPE.Placebo)
inexact.power.std.RMSPE.Placebo = apply(inexact.powers.RMSPE.Placebo, 1, sd)

inexact.power.mean.over.time.Placebo = apply(inexact.powers.over.time.Placebo, c(1,2), mean)
inexact.power.std.over.time.Placebo = apply(inexact.powers.over.time.Placebo, c(1,2), sd)

### Get randomized placebo results
random.power.mean.RMSPE.Placebo = rowMeans(random.powers.RMSPE.Placebo)
random.power.std.RMSPE.Placebo = apply(random.powers.RMSPE.Placebo, 1, sd)

random.power.mean.over.time.Placebo = apply(random.powers.over.time.Placebo, c(1,2), mean)
random.power.std.over.time.Placebo = apply(random.powers.over.time.Placebo, c(1,2), sd)
######################################################################################
######################################################################################
######################################################################################
### Comparison ###

library(reshape2)
df1 <- data.frame(#power.mean.RMSPE.Placebo,
                  random.power.mean.RMSPE.Placebo,
                  inexact.power.mean.RMSPE.Placebo,
                  power.mean.RMSPE.LOO, 
                  power.mean.RMSPE.LTO,
                  power.minusc.mean.RMSPE.LTO,
                  taus)
df2 <- melt(df1, id.vars='taus')


sds = melt(data.frame(#power.std.RMSPE.Placebo,
                      random.power.mean.RMSPE.Placebo,
                      inexact.power.std.RMSPE.Placebo,
                      power.mean.RMSPE.LOO,
                      power.std.RMSPE.LTO,
                      power.minusc.std.RMSPE.LTO,
                      taus), id.vars='taus')
df2$sd = sds$value
df2

#################################################
############ Set Colors #################
#################################################

# colors via hcl.colors
methodColors <-
  setNames( c('#a6cee3','#1f78b4','#00366C',"#D4F3A3","#65BB8C","#0D7434")
            , c('power.mean.RMSPE.Placebo',
                'inexact.power.mean.RMSPE.Placebo',
                'random.power.mean.RMSPE.Placebo',
                'power.mean.RMSPE.LOO',
                'power.mean.RMSPE.LTO',
                'power.minusc.mean.RMSPE.LTO')  )

methodNames <-   setNames( c('Naive Placebo','Appr. Placebo', 'Rand. Placebo', 'LOO', 'Naive LTO','Powered LTO')
                           , c('power.mean.RMSPE.Placebo',
                               'inexact.power.mean.RMSPE.Placebo',
                               'random.power.mean.RMSPE.Placebo',
                               'power.mean.RMSPE.LOO',
                               'power.mean.RMSPE.LTO',
                               'power.minusc.mean.RMSPE.LTO')  )


#################################################
############ Section 4 Plots #################
#################################################

library(ggplot2)
ggplot(df2, aes(x=as.factor(taus), y=value, fill=variable)) +
  geom_bar(stat='identity', color = 'black', position=position_dodge()) +
  #geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width = 0.2,
  #              position=position_dodge(.9)) + 
  scale_fill_manual(values = methodColors,name = "Method", labels=methodNames) +
  theme_bw() +
  ggtitle("Prop. 99 Data: RMSPE, N = 30, alpha = 0.02") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Effect Sizes (standard deviations)") + 
  ylab("Power") +
  scale_x_discrete(labels= c(-3,-2,-1,0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.9)) +
  geom_hline(yintercept = 0.02, linetype = 2, color = "red")


#plotname = sprintf("../Figures/smoking_RMSPE_N%d_alpha%.2f.pdf", size.resampled.dataset,alpha )
#plotname = "smoking_rmspe_n30_alpha0.02.png"
dev.print(png, plotname, width = 500, height = 400)

#################################################
############ Introduction Plots #################
#################################################

methodNames <-   setNames( c('Naive Placebo','Appr. Placebo','Rand. Placebo', 'LTO','Powered LTO')
                           , c('power.mean.RMSPE.Placebo',
                               'inexact.power.mean.RMSPE.Placebo',
                               'random.power.mean.RMSPE.Placebo',
                               'power.mean.RMSPE.LTO',
                               'power.minusc.mean.RMSPE.LTO')  )

df1 <- data.frame(#power.mean.RMSPE.Placebo, 
  inexact.power.mean.RMSPE.Placebo,
  power.mean.RMSPE.LTO,
  #power.minusc.mean.RMSPE.LTO,
  taus)
df2 <- melt(df1, id.vars='taus')


sds = melt(data.frame(#power.std.RMSPE.Placebo,
  inexact.power.std.RMSPE.Placebo,
  power.std.RMSPE.LTO,
  #power.minusc.std.RMSPE.LTO,
  taus), id.vars='taus')
df2$sd = sds$value
df2

df2$variable <- as.factor(df2$variable)


ggplot(df2, aes(x=as.factor(taus), y=value, fill=variable )) +
  geom_bar(stat='identity', color = 'black', position=position_dodge()) +
  #geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width = 0.2,
  #              position=position_dodge(.9)) + 
  #scale_fill_discrete(labels=c('Exact Placebo', 'Appr. Placebo', 'Naive LTO', 'Powered LTO')) + +
  scale_fill_manual(values = methodColors,name = "Method", labels=methodNames) +
  theme_bw() +
  ggtitle("Prop. 99 Data: RMSPE, N = 30, alpha = 0.02") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Effect Sizes (standard deviations)") + 
  ylab("Power") +
  scale_x_discrete(labels= c(-3,-2,-1,0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.9))+
  geom_hline(yintercept = 0.02, linetype = 2, color = "red")

#plotname = sprintf("../Figures/smoking_RMSPE_N%d_alpha%.2f.pdf", size.resampled.dataset,alpha )
plotname = "smoking_rmspe_n30_alpha0.02_intro.png"
dev.print(png, plotname, width = 500, height = 400)

#################################################
### Power at a fixed time
#################################################

###############################################################
############# Plotting Colors #################################
###############################################################
time = 31
df1 <- data.frame(power.mean.over.time.Placebo[,time], 
                  inexact.power.mean.over.time.Placebo[,time],
                  power.mean.over.time.LTO[,time],
                  power.minusc.mean.over.time.LTO[,time],
  taus)


df2 <- melt(df1, id.vars='taus')


sds = melt(data.frame(power.std.over.time.Placebo[,time],
                      inexact.power.std.over.time.Placebo[,time],
                      power.std.over.time.LTO[,time],
                      power.minusc.std.over.time.LTO[,time],
  taus), id.vars='taus')
df2$sd = sds$value
df2

methodColors <-
  setNames( c('#a6cee3','#1f78b4','#b2df8a','#33a02c')
            , unique(df2$variable)  )

methodNames <-   setNames( c('Naive Placebo','Appr. Placebo', 'Naive LTO','Powered LTO')
                           , unique(df2$variable)  )

###############################################################
############# Plotting #################################
###############################################################

time = 31
df1 <- data.frame(#power.mean.over.time.Placebo[,time], 
                  inexact.power.mean.over.time.Placebo[,time],
                  power.mean.over.time.LTO[,time],
                  power.minusc.mean.over.time.LTO[,time],
                  taus)

df2 <- melt(df1, id.vars='taus')


sds = melt(data.frame(#power.std.over.time.Placebo[,time],
                      inexact.power.std.over.time.Placebo[,time],
                      power.std.over.time.LTO[,time],
                      power.minusc.std.over.time.LTO[,time],
                      taus), id.vars='taus')
df2$sd = sds$value
df2

library(ggplot2)
ggplot(df2, aes(x=as.factor(taus), y=value, fill=variable)) +
  geom_bar(stat='identity', color = 'black', position=position_dodge()) +
  #geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width = 0.2,
  #              position=position_dodge(.9)) + 
  scale_fill_manual(values = methodColors, name = "Method", labels=methodNames) +
  theme_bw() +
  ggtitle("Prop. 99 Data: Year 2000 difference, N = 30, alpha = 0.02") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Effect Sizes (standard deviations)") + 
  ylab("Power") +
  scale_x_discrete(labels= c(-3,-2,-1,0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) 

#plotname = sprintf("../Figures/smoking_fixed_time_N%d_alpha%.2f.pdf", size.resampled.dataset,alpha )
plotname = "smoking_fixedtime_n30_alpha0.02.png"
dev.print(png, plotname, width = 500, height = 400)


### Power over time ###

### Higher Tau ###
time = 27:31
idx = 4
df1 <- data.frame(#power.mean.over.time.Placebo[idx,time], 
                  inexact.power.mean.over.time.Placebo[idx,time],
                  power.mean.over.time.LTO[idx,time],
                  power.minusc.mean.over.time.LTO[idx,time],
                  time)
df2 <- melt(df1, id.vars='time')


sds = melt(data.frame(#power.std.over.time.Placebo[idx,time],
                      inexact.power.std.over.time.Placebo[idx,time],
                      power.std.over.time.LTO[idx,time],
                      power.minusc.std.over.time.LTO[idx,time],
                      time), id.vars='time')
df2$sd = sds$value
df2

ggplot(df2, aes(x=as.factor(time), y=value, fill=variable)) +
  geom_bar(stat='identity', color = 'black', position=position_dodge()) +
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width = 0.2,
                position=position_dodge(.9))+ 
  scale_fill_discrete(name = "Method",labels=c('Placebo', 'LTO', 'Powered LTO')) +
  ggtitle(sprintf("Smoking Data: Effect Size: %.1f, N = 15, alpha = 0.02",taus[idx])) +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Year") + 
  ylab("Power")+
  scale_x_discrete(labels= 1969+time) 
  #coord_cartesian(ylim=c(0,0.04))

plotname = sprintf("../Figures/smoking_overtime_N%d_alpha%.2f_tau%.2f.pdf", size.resampled.dataset,alpha,taus[idx] )
dev.print(pdf,plotname)
