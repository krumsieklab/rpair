rm(list = ls())
library(glmnet)
library(magrittr)

df = sgi::qmdiab_clin
X = sgi::qmdiab_plasma %>% as.matrix() %>% scale

y = df$Diabetes %>% as.numeric()

set.seed(42)
cv = cv.glmnet(X[1:200,], y[1:200], nfolds = 5, keep = TRUE)
kc.cv = kc.cv.glmnet(X[1:200,], y[1:200], nfolds = 5, keep = TRUE)

plot(cv)
abline(v = log(cv$lambda.min), col = "blue" )
abline(v = log(cv$lambda.1se), col = "green" ) # default
abline( h = cv$cvm[cv$lambda==cv$lambda.1se], col = "green")

# default option for glmnet cv is lmabda.1se
coef(cv) %>% as.vector() %>% {.[.!=0]}
# lambda.min coefficient
coef(cv, s = "lambda.min") %>% as.vector() %>% {.[.!=0]}


library(survival)
# pair example 
xh = X[1:100,1:10]
S = Surv(seq(100), seq(100)*0+1 )
foldids = sample(100)

# train-test set 1
xh_tr = xh[foldids!=1,]
S_tr = S[foldids!=1,]

# continuous outcome
y = rnorm(100)

# train-test set 1
y_tr = y[foldids!=1]


# pairs 
pi = cbind(sample(100, 500, replace = T), sample(100, 500, replace = T))
pi_tr =  pi[ !(foldids[pi[,1]] == 1 | foldids[pi[,2]] == 1),]