Tutorial 4: supported outcome types
================

This tutorial showcases each of the different types of outcome / response variables supported by the rpair. These are: 
survival outcome (2 columns), survival outcome (3 columns), numeric outcome, ordinal outcome, factor outcome, and 
comparable pairs. Each example below creates a simulated dataset of the corresponding type.


``` r
library(rpair)
library(magrittr)
library(survival)
```

## Survival outcome (2 columns)
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
# generate survival times
S = Surv(exp(scale(X %*% b)), sample(c(F,T),ks, T) )
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
# @mubu - generate random data with 3-column survival outcome
```

## Numeric outcome
```r
# @mubu - generate random data with continuous outcome
```

## Ordinal outcome
```r
# @mubu - generate random data with ordinal outcome
```

## Comparable pairs
Users can also provide comparable pairs as direct input to the cv_rpair (and rpair) function(s). The following
example takes the 2-column survival outcome produced in the first example and uses the internal rpair function to
generate comparable pairs.

IMPORTANT NOTE: The function cv_rpair does not support the use of comparable pairs without user-provided folds. This is because important information for generating stratified folds (for example, in the case of comparable pairs generated from survival data) can be lost when only providing the pairs.
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

# generate survival times
S = Surv(exp(scale(X %*% b)), sample(c(F,T),ks, T) )
colnames(S) = c("time", "status")

# generate comparable pairs from 2-column survival data
cp = rpair:::y_to_pairs(S)

# get fold ids
fids = rpair:::get_stratified_folds(S, nfolds=5)

# fit model
surv2_fit = cv_rpair(X, cp, foldid=fids$fold)

# get log lambda plot and beta coefficients
# @ KC - add plot and coef() results
```
