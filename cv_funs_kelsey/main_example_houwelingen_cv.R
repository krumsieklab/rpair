# implementation of Houwelingen type CV loss calculation

rm(list=ls())
library(survival)
library(magrittr)
library(glmnet)

source("cv_funs_kelsey/cv_rpair.R")
source("rpair_gloss.R")
source("rpair_loss_eval_funs.R")
source("cv_funs_kelsey/ssvm_utils.R")
source("utils_pairs.R")
source("utils_repair.R")
source("load_fortran.R")

# generate example data  --------------------------------------------------

# number of samples 
ks = 100 

# number of features
n1 = 20 # n informative features
n2 = 40 # n noisy features

set.seed(42)

# coefficients 
b = c( exp(runif(n1)), seq(n2)*0) %>% round(2)

x = matrix( rnorm((n1+n2)*100), nrow = ks ) %>% scale
# survival outcome
S = Surv(exp( scale(x %*% b)), sample(c(F,T),ks, T) )
colnames(S) = c("time", "status")

# fold ids
fids = sample(4, ks, T)
#----------------------------------------------------------------

xtr = x
Str = S

alpha = 1
# the maximum number of variables ever to be nonzero 
pmx = min(ncol(x), sum(S[,2]))
if(alpha < 0.5 ) pmx = ncol(x)+1
loss_type = "exp" # "log"

kpcv = cv.rpair(xtr, Str, foldid=fids, nlambda=25, type.measure = "deviance", alpha = alpha, 
                pmax = pmx, alignment = "fraction", keep = T, loss_type = loss_type)


# all pairs
pi_all = y_to_pairs(Str)

# there are two loss functions 
exp_loss <- function(x) exp(x)
log_loss <- function(x) log(1+exp(x))
ffs <- list(exp = exp_loss, log = log_loss)

# losses 
lo = 
sapply(seq(unique(fids)), function(i){
  # fold i 
  # training set samples
  j = which(fids != i)
  
  # all_pairs - training_set_pairs
  pir =  pi_all[!( (pi_all[,1] %in% j) & (pi_all[,2] %in% j) ),]
  print(nrow(pir))

  # loss margin z i.e. loss_function(z)
  Z = kpcv$fit.preval[[i]][pir[,1],] -  kpcv$fit.preval[[i]][pir[,2],]
  # Z is a matrix of size:n_valid_pairs X n_lambda
  
  # pick loss function based on the model loss type 
  ff = ffs[[loss_type]] # ffs is the list of functions 
  
  # for each lambda calculate deviance for fold i 
  # colMeans(exp_loss(Z),na.rm = T)
  colMeans( ff(Z), na.rm = T )
})

# that is the estimated CV deviance i.e. mean over folds 
deviance =  rowMeans(lo,na.rm = T) 

plot(log(deviance))
abline(v = which.min(deviance), col = "red", lty = 2)

# coefficients
kpcv$rpair.fit$beta[,which.min(deviance)] %>% plot
