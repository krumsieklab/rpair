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
library(magrittr)
library(ggplot2)
library(tidyr)
library(dplyr)

rm(list = ls())

# where latest code is 
#setwd("~/Documents/R/p_repair/")
setwd("~/rpair/")

# utility functions for pair 
source("utils_pairs.R")

# utility functions for repair model fit 
source("utils_repair.R")

# load backend fortran workhorse
source("load_fortran.R")

# core model fitting for exp and log 
source("rpair_gloss.R")
source("repair_gloss.R")

# evaluation functions for rpair objects
source("rpair_loss_eval_funs.R")

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
k_fite = rpair_gloss(x, cp, standardize = F, pmax = 50, loss_type = "exp")
m_fite = repair_gloss(x, cp, standardize = F, pmax = 50, loss_type = "exp")
# logistics loss
k_fitl = rpair_gloss(x, cp, standardize = F, pmax = 50, loss_type = "log")
m_fitl = repair_gloss(x, cp, standardize = F, pmax = 50, loss_type = "log")

k_fite$lambda
k_fite$beta[,10]
k_fite$dev.ratio %>%  plot

plot(k_fite)


# test predict function
g_fite <- k_fite
class(g_fite) <- "glmnet"
g_fite$a0 <- rep(0, g_fite$dim[2])
g_fite$offset <- FALSE

k_tmp <- predict(k_fite, type = "coefficients")
View(k_tmp %>% as.matrix())
g_tmp <- predict(g_fite, type = "coefficients")
View(g_tmp %>% as.matrix())

k_tmp <- predict(k_fite, type = "nonzero")
View(k_tmp)
g_tmp <- predict(g_fite, type = "nonzero")
View(g_tmp)

k_tmp <- predict(k_fite, x, type="link")
View(k_tmp)
g_tmp <- predict(g_fite, x, type="link")
View(g_tmp)
identical(k_tmp, g_tmp)


# test coef function
k_tmp <- coef(k_fite)
View(k_tmp)
g_tmp <- coef(g_fite)


# test plot function
k_plot <- plot(k_fite, xvar="norm")
plot(k_plot)
g_plot <- plot(g_fite, xvar="norm")
k_plot <- plot(k_fite, xvar="lambda", label=TRUE)
plot(k_plot)
g_plot <- plot(g_fite, xvar="lambda")
k_plot <- plot(k_fite, xvar="dev", legend=TRUE)
plot(k_plot)
g_plot <- plot(g_fite, xvar="dev")






