Tutorial 4: supported outcome types
================

This tutorial showcases each of the different types of outcome / response variables supported by the rpair. These are: 
survival outcome (2 columns), survival outcome (3 columns), numeric outcome, ordinal outcome, factor outcome, and 
comparable pairs. For this example, we first create a single simulated dataset and a vector of coefficients. We then
use these to generate the specific outcomes in each of their corresponding sections.


``` r
library(rpair)
library(magrittr)
library(survival)
```

## Generate Random Data
```r
set.seed(73)

# generate survival data
ks = 100

# number of features
n1 = 20 # n informative features
n2 = 40 # n noisy features

# coefficients
b = c(exp(runif(n1)), seq(n2)*0) %>% round(2)

# generate data matrix
X = matrix(rnorm((n1+n2)*100), nrow = ks) %>% scale
```

## Survival outcome (2 columns)
```r
# generate two-column survival object
S = Surv(exp(scale(X %*% b))*10, sample(c(F,T),ks, T) )
colnames(S) = c("time", "status")

surv2_fit = cv_rpair(X, S)

plot(surv2_fit)

coef(surv2_fit)[1:10]

```
<img src="imgs/surv2_plot.png" width="665" height="455" />

[1]  0.000000000  0.000000000 -0.416909156 -0.028875330 -0.027832046 -0.003926241  0.000000000  0.000000000  0.000000000   
[10] -0.162762181


## Survival outcome (3 columns e.g. count type Surv object)
```r
# generate three-column surival object
S = Surv(rep(0, ks), exp(scale(X %*% b))*10, sample(c(F,T),ks, T) )
#colnames(S) = c("start", "stop", "status")

surv3_fit = cv_rpair(X, S)

plot(surv3_fit)

coef(surv3_fit)

pred(surv3_fit, X)
```

## Numeric outcome
```r
# numeric outcome here is same as time and stop for two survival outcomes
y = exp(scale(X %*% b))*10

num_fit = cv_rpair(X, y)

plot(num_fit)

coef(num_fit)

pred(num_fit, X)
```

## Ordinal outcome
```r
# create ordinal outcome 
y = scale(X %*% b)
y = (y - min(y) + 10)^2
# number of levels for ordinal outcome
l = 4
# ordinal levels
y = cut(y, breaks = seq(min(y), max(y), length.out = l+1), include.lowest = T) %>% as.numeric

ord_fit = cv_rpair(X, y)

plot(ord_fit)

coef(ord_fit)

pred(ord_fit, X)
```

## Comparable pairs
Users can also provide comparable pairs as direct input to the cv_rpair (and rpair) function(s). The following
example takes the 2-column survival outcome produced in the first example and uses the internal rpair function to
generate comparable pairs.

IMPORTANT NOTE: The function cv_rpair does not support the use of comparable pairs without user-provided folds. This is because important information for generating stratified folds (for example, in the case of comparable pairs generated from survival data) can be lost when providing only the pairs.

```r
# generate survival times
S = Surv(exp(scale(X %*% b)), sample(c(F,T),ks, T) )
colnames(S) = c("time", "status")

# generate comparable pairs from 2-column survival data
cp = rpair:::y_to_pairs(S)

# get fold ids
fids = rpair:::get_stratified_folds(S, nfolds=5)

# fit model
cp_fit = cv_rpair(X, cp, foldid=fids$fold)

# get log lambda plot and beta coefficients
plot(cp_fit)

coef(cp_fit)

pred(cp_fit, X)
```
