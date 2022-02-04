rm(list=ls())
library(survival)
library(magrittr)
library(glmnet)

setwd("~/Downloads/rstudio-export_kelsey/")
source("fishnet.R")
source("lognet.R")
source("ssvm_utils.R")

set.seed(42)
load("rppa_nods")
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

# benchhorse function
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

# debug the two manually implemented cv functions
pcv = cv.srpoisson(xtr, Str, fids,k = 25, type.measure = "deviance", alpha = alpha, pmx = pmx)
# fishnet line 20
# Warning message:
#  from glmnet Fortran code (error code -10018); Number of nonzero coefficients along the path exceeds pmax=56 
#  at 18th lambda value; solutions for larger lambdas returned 


bcv = cv.srbinomial(xtr, Str, fids, k = 25, type.measure = "auc", alpha = alpha, pmx = pmx)
# lognet line 30
#Warning messages:
#  1: In lognet(xd, is.sparse, ix, jx, y, weights, offset, alpha, nobs,  :
#     one multinomial or binomial class has fewer than 8  observations; dangerous ground
#  2: from glmnet Fortran code (error code -10019); Number of nonzero coefficients along the path exceeds pmax=56 
#     at 19th lambda value; solutions for larger lambdas returned 

