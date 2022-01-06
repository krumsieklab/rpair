# ff(
#   x, # data 
#   y, # outcome 
#   family, # we have `loss` instead,
#   weights=NULL, # x,
#   offset=NULL, # x~,
#   alpha=1.0, # yes
#   nlambda=100, # yes
#   lambda.min.ratio=ifelse(nobs<nvars,1e-2,1e-4), # yes
#   lambda=NULL, # yes,
#   standardize=TRUE, #yes?
#   intercept=TRUE, # x
#   thresh=1e-7, # yes 
#   dfmax=nvars+1, # yes?
#   pmax=min(dfmax*2+20,nvars), # yes 
#   exclude=NULL, # x
#   penalty.factor=rep(1,nvars), # yes 
#   lower.limits=-Inf, # ?? GGM 
#   upper.limits=Inf, # ??
#   maxit=100000, # yes 
#   type.gaussian=ifelse(nvars<500,"covariance","naive"), # x 
#   type.logistic=c("Newton","modified.Newton"), # yes 
#   standardize.response=FALSE, # x
#   type.multinomial=c("ungrouped","grouped"), # x
#   relax=FALSE, # x
#   trace.it=0, # yes 
#   ...){
#   


# log -exp losses are working fine 
# with a caveat that log-loss does not accept defined lambda 
# sequence, it always insists to calculate the seq each time 
rm(list = ls())

# where latest code is 
setwd("~/Documents/R/p_repair/")

# utility functions for pair 
source("utils_pairs.R")

# utility functions for repair model fit 
source("utils_repair.R")

# load backend fortran workhorse
source("load_fortran.R")

# core model fitting for exp and log 
source("repair_gloss.R")

# convert to pairs 
fp<- function(S){
  time <- S[,1]
  status <- S[,2]
  N = length(time)
  # for tied times
  time[status == 0] = time[status == 0]+1e-4
  dtimes <- time  
  dtimes[status == 0] = Inf
  which(outer(time, dtimes, ">"), arr.ind = T)
}

# generate some random data 
set.seed(41)
x = matrix(rnorm(40000),ncol = 200 )
S = cbind(sample(nrow(x)), rbinom(nrow(x),1,prob = 0.7))

cp = fp(S)
# exponential loss
fite = repair_gloss(x, cp, standardize = F, pmax = 50, loss_type = "exp")
# logistics loss
fitl = repair_gloss(x, cp, standardize = F, pmax = 50, loss_type = "log")

fite$lambda
fite$beta[,10]
fite$dev.ratio %>%  plot

plot(fite)







