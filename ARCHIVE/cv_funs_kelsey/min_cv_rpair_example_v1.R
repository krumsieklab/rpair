rm(list=ls())
library(survival)
library(magrittr)
library(glmnet)

PATH_TO_INPUT = "~/Desktop/rppa_nods"

setwd("~/rpair/")
source("cv_funs_kelsey/cv_rpair.R")
source("rpair_gloss.R")
source("rpair_loss_eval_funs.R")
source("cv_funs_kelsey/ssvm_utils.R")
source("utils_pairs.R")
source("utils_repair.R")
source("load_fortran.R")

set.seed(42)
load(PATH_TO_INPUT)
# remove missing values
dfnods = dfnods %>% lapply( na.omit)

# train test splits and inner cvs
dfnods<- dfnods %>% lapply( function(nod, nboo=100) 
  list( data = nod, spids = replicate(nboo, {
    ids = fcvid(nod$S, 5) < 5
    ids[ids] = fcvid(nod$S[ids,], 5)
    ids
  })))

poc <- dfnods$MESO

# set up the data
x = as.matrix( scale( poc$data$rppa ) )
S = poc$data$S

# CV folds for each tr
spids <- poc$spids
rm(poc)

# take one of the folds to test
trinds = spids[,1]
alpha=1
# weird glmnet error for surv.time == 0
S[,1] = S[,1] + 1

# training data 1,2,3,4,5:train data
xtr = x[trinds>0, ]
Str = S[trinds>0, ]

# test data 0:test data
Sts = S[trinds==0, ]
xts = x[trinds==0, ]

# cv ids
fids = trinds[trinds>0] 

# the maximum number of variables ever to be nonzero 
pmx = min(ncol(x), sum(S[,2]))
if(alpha < 0.5 ) pmx = ncol(x)+1

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
                pmax = pmx, alignment = "lambda", nlambda = 25, keep = T, loss_type="log")
# all fit.preval dfs are equal entirely 0
# gives thhese warning for every fold:
#  - from glmnet Fortran code (error code -10001); Number of nonzero coefficients along the path exceeds pmax=56 at 
#     1th lambda value; solutions for larger lambdas returned
#  - In getcoef.glmnetfit(fit, nvars, nx, vnames) :
#       an empty model has been returned; probably a convergence issue

