library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

setwd("~/git/censoredShiftedWald") # update this to your local repo

# load datasets

# Data Set 1 - Faulkenberry, Vick, & Bowman (2018)
data1 = data.frame()
filestem = "https://raw.githubusercontent.com/tomfaulkenberry/physNumComparisonTask/master/results/data/subject_"
for (n in 101:123){
  file = paste(filestem, n, ".csv", sep="")
  temp = read.csv(file)
  temp$subject_nr = n - 100
  data1 = rbind(data1,temp)
}  

# Data Set 2 -- Bowman & Faulkenberry (2020)
data2 = data.frame()
filestem = "https://raw.githubusercontent.com/Kbow27/Thesis/master/subject-"
for (n in 101:153){
  file = paste(filestem, n, ".csv", sep="")
  temp = read.csv(file)
  temp$subject_nr = (n - 100) + 23
  data2 = rbind(data2, temp)
}

# Data Set 3 -- unpublished data set
data3 = data.frame()
filestem = "https://raw.githubusercontent.com/tomfaulkenberry/physNumComparisonTask/master/results/dataset3/subject-"
for (n in 1:35){
  file = paste(filestem, n, ".csv", sep="")
  temp = read.csv(file)
  temp$subject_nr = n + 23 + 53
  data3 = rbind(data3,temp)
}

raw = rbind(data1,data2,data3)

dat1 = raw %>%
  subset(congruity == "congruent")

dat2 = raw %>%
  subset(congruity == "incongruent")

nSub1 = length(unique(dat1$subject_nr))
nTrials1 = length(dat1$response_time[dat1$subject_nr==1])
nSub2 = length(unique(dat2$subject_nr))
nTrials2 = length(dat2$response_time[dat2$subject_nr==1])

# set up model inputs
# key: 1 = congruent, 2 = incongruent

rts1 = matrix(0, nrow=nSub1, ncol=nTrials1)
err1 = matrix(0, nrow=nSub1, ncol=nTrials1)

for (i in 1:nSub1){
  rts1[i,] = dat1$response_time[dat1$subject_nr==i]/1000
  err1[i,] = 1 - dat1$correct[dat1$subject_nr==i]
}

# handle outliers by censoring
shortThresh = 0.2
longThresh = 2
rts1[rts1 < shortThresh] <- shortThresh
rts1[rts1 > longThresh] <- longThresh

rts2 = matrix(0, nrow=nSub2, ncol=nTrials2)
err2 = matrix(0, nrow=nSub2, ncol=nTrials2)

for (i in 1:nSub2){
  rts2[i,] = dat2$response_time[dat2$subject_nr==i]/1000
  err2[i,] = 1 - dat2$correct[dat2$subject_nr==i]
}

# handle outliers by censoring
shortThresh = 0.2
longThresh = 2
rts2[rts2 < shortThresh] <- shortThresh
rts2[rts2 > longThresh] <- longThresh


# fit censored shifted Wald

# fit model using HMC
fit1 = stan(file = "swald.stan", 
            data = list(rt = rts1,
                        D = err1,
                        nrt = nTrials1, 
                        ns = nSub1)
)

fit2 = stan(file = "swald.stan", 
            data = list(rt=rts2,
                        D = err2,
                        nrt=nTrials2, 
                        ns=nSub2)
)

# extract estimates
estimates = matrix(0, nrow=nSub1, ncol=6)

for (i in 1:nSub1){
  estimates[i,1] = summary(fit1, pars=c(sprintf("G[%s]",i)))$summary[,"mean"]
  estimates[i,2] = summary(fit2, pars=c(sprintf("G[%s]",i)))$summary[,"mean"]
  estimates[i,3] = summary(fit1, pars=c(sprintf("A[%s]",i)))$summary[,"mean"]
  estimates[i,4] = summary(fit2, pars=c(sprintf("A[%s]",i)))$summary[,"mean"]
  estimates[i,5] = summary(fit1, pars=c(sprintf("H[%s]",i)))$summary[,"mean"]
  estimates[i,6] = summary(fit2, pars=c(sprintf("H[%s]",i)))$summary[,"mean"]
}

estimates = as.data.frame(estimates)
names(estimates) <- c("drift1","drift2", "thresh1", "thresh2", "ndt1","ndt2")
write.csv(estimates, "estimates.csv", row.names=F)


