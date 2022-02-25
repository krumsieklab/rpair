# Houwelingen loss example --------------------------------------------------

rm(list=ls())
library(survival)
library(magrittr)
library(rpair)

# generate example data  --------------------------------------------------

# number of samples
ks = 100

# number of features
n1 = 10 # n informative features
n2 = 300 # n noisy features

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

kpcv5 = cv_rpair(xtr, Str, foldid=fids, nlambda=25, type.measure = "cindex", alpha = alpha,
                pmax = pmx, alignment = "fraction", keep = T, loss_type = loss_type)

kpcv6 = cv_rpair(xtr, Str, foldid=fids, nlambda=25, type.measure = "cindex", alpha = alpha,
                 pmax = pmx, alignment = "fraction", keep = T, loss_type = "log")

plot(kpcv5)
plot(kpcv6)

