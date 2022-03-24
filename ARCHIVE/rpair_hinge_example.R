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
# squared hinge loss
sfit = rpair_hinge(x, cp, standardize = F, pmax = 50, loss_type = "sqh")
# huberized squared hinge loss
hfit = rpair_hinge(x, cp, standardize = F, pmax = 50, loss_type = "huh")

# predict.rpair
plot(hfit, xvar="norm")
splot <- plot(sfit, xvar="lambda")
plot(splot)

# coef.rpair
hcoef <- coef(hfit)
View(hcoef %>% as.matrix())
scoef <- coef(sfit)
View(scoef %>% as.matrix())

# plot.rpair
pcoef <- predict(sfit, type = "coefficients")
View(pcoef %>% as.matrix())
pnz <- predict(hfit, type = "nonzero")
View(pnz)
pred <- predict(sfit, x, type="link")
View(pred)

