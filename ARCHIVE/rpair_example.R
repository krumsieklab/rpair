# log -exp losses are working fine
# with a caveat that log-loss does not accept defined lambda
# sequence, it always insists to calculate the seq each time
library(rpair)

rm(list = ls())

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
efit = rpair_gloss(x, cp, standardize = F, pmax = 50, loss_type = "exp")
# logistics loss
lfit = rpair_gloss(x, cp, standardize = F, pmax = 50, loss_type = "log")

# predict.rpair
plot(lfit, xvar="norm")
lplot <- plot(efit, xvar="lambda")
plot(lplot)
dplot <- plot(lfit, xvar="dev", legend=TRUE)
plot(dplot)

# coef.rpair
ecoef <- coef(efit)
View(ecoef %>% as.matrix())
lcoef <- coef(lfit)
View(lcoef %>% as.matrix())

# plot.rpair
pcoef <- predict(lfit, type = "coefficients")
View(pcoef %>% as.matrix())
pnz <- predict(efit, type = "nonzero")
View(pnz)
pred <- predict(lfit, x, type="link")
View(pred)

