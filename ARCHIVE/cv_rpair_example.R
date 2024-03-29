# Example 2 (from Houwelingen loss example script)

rm(list=ls())
library(survival)
library(rpair)


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

# predict.cv_rpair
pred1 = predict(kcv1, newx = xtr, s="lambda.min")
pred2 = predict(kcv2, newx = xtr, s="lambda.min")
pred3 = predict(kcv3, newx = xtr, s="lambda.1se")
pred4 = predict(kcv4, newx = xtr, s="lambda.1se")

# coef.cv_rpair
coef1 = coef(kcv1, s="lambda.1se")
coef2 = coef(kcv2, s="lambda.1se")
coef3 = coef(kcv3, s="lambda.min")
coef4 = coef(kcv4, s="lambda.min")

# plot.cv_rpair
plot(kcv1, ggplot=F)
plot(kcv2)
plot(kcv3, ggplot=F)
plot(kcv4)


