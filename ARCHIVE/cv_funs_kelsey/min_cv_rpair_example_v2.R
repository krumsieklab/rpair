# Example 2 (from Houwelingen loss example script)

rm(list=ls())
library(survival)
library(magrittr)
library(glmnet)

setwd("~/rpair/")

source("cv_funs_kelsey/cv_rpair.R")
source("rpair_gloss.R")
source("rpair_loss_eval_funs.R")
source("cv_funs_kelsey/ssvm_utils.R")
source("utils_pairs.R")
source("utils_repair.R")
source("load_fortran.R")
source("cv_funs_kelsey/cv_utils.R")

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

# exp loss
kcv1 = cv_rpair(xtr, Str, foldid=fids, nlambda=25, type.measure = "deviance", alpha = alpha, 
                pmax = pmx, alignment = "fraction", keep = T, loss_type="exp")

# log loss
kcv2 = cv_rpair(xtr, Str, foldid=fids, type.measure = "deviance", alpha = alpha, 
                pmax = pmx, alignment = "fraction", nlambda = 25, keep = T, loss_type="log")

# fold-ids missing
kcv3 = cv_rpair(xtr, Str, nfolds = 5, type.measure = "deviance", alpha = alpha, 
                pmax = pmx, alignment = "fraction", nlambda = 25, keep = T, loss_type="log")

# alignment = lambda - this does not appear to be working correctly
kcv4 = cv_rpair(xtr, Str, foldid = fids, type.measure = "deviance", alpha = alpha, 
                pmax = pmx, alignment = "lambda", nlambda = 25, keep = T, loss_type="exp")
# all fit.preval dfs are equal entirely 0
# gives this warning for every fold:
#  from glmnet Fortran code (error code -10001); Number of nonzero coefficients along the path exceeds pmax=48 at 
#    1th lambda value; solutions for larger lambdas returned 


#----------------------------------------------------------------
kcoef <- coef(kcv4)
kpred <- predict(kcv4, xtr)
plot(kcv4, ggplot=T)

