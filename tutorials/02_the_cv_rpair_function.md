Tutorial 2: the cv_rpair function
================

This tutorial demonstrates how to use the cv_rpair function and its utility functions. The cv_rpair function performs
k-fold cross-validation of the rpair function. The cv_rpair functions supports two kinds of input: survival and continuous. Like rpair, cv_rpair has three corresponding utility functions: predict, coef, and plot. 
An important parameter featured in this tutorial is type.measure. The type.measure parameter allows  the user to select which loss function to use for cross-validation - "deviance" (default) or "cindex". Selecting "deviance" uses the selected loss function itself (e.g. exponential loss) while selecting "cindex" uses the concordance measure.

``` r
library(rpair)
library(magrittr)
library(survival)
```

## Generate Random Survival Dataset
```r
# number of samples
ks = 100

# number of features
n1 = 20 # n informative features
n2 = 40 # n noisy features

set.seed(42)

# coefficients
b = c( exp(runif(n1)), seq(n2)*0) %>% round(2)

x = matrix( rnorm((n1+n2)*100), nrow = ks ) %>% scale
S = Surv(exp( scale(x %*% b)), sample(c(F,T),ks, T) )
colnames(S) = c("time", "status")

fids = sample(4, ks, T)

alpha = 1
# the maximum number of variables ever to be nonzero
pmx = min(ncol(x), sum(S[,2]))
if(alpha < 0.5 ) pmx = ncol(x)+1


```

```r
# first 4 rows and columns of the data matrix
x[1:4,1:4]
```
               [,1]       [,2]        [,3]        [,4]
    [1,]  1.2613010  0.1289885 -1.29395115 -0.08864345
    [2,]  2.2108125  0.2711029  2.08074800 -0.84310551
    [3,] -1.3439050 -0.3622906  0.07022594 -1.14517033
    [4,] -0.2703135 -0.3823347 -0.12974566 -0.70252506 

```r
# first five rows of survival data
colnames(S) <- c("time", "status")
S[1:5,]
```
    [1] 0.2499683  7.2891822  2.3499586  1.6015316+ 1.2846922 

## type.measure = "deviance"
```r
efit = cv_rpair(x, S, foldid=fids, nlambda=25, type.measure = "deviance", alpha = alpha,
                pmax = pmx, alignment = "fraction", keep = T, loss_type="exp")
plot(efit)
```

<img src="imgs/cv_efit_plot.png" width="665" height="455" />

```r
epred = predict(efit, newx = x, s="lambda.min")
epred[1:5]
```
    [1]  2.0677687 -2.9669107 -1.0780859 -1.3149500 -0.5849274
    
```r
ec = coef(efit, s="lambda.1se")
ec[1:5]
```
    [1] -0.33103864 -0.09459406  0.00000000  0.00000000 -0.03626532

```r
hfit = cv_rpair(x, S, loss_type="huh", foldid=fids, nlambda=25, type.measure = "deviance",
                pmax = pmx, keep = T, alignment="lambda")
plot(hfit, ggplot=F)
```

<img src="imgs/cv_hfit_plot.png" width="665" height="455" />

```r
hpred = predict(hfit, x, s="lambda.1se")
hpred[1:5]
```
    [1]  0.4348660 -0.8906891 -0.3107375 -0.4157110 -0.1552527

```r
hc = coef(hfit, s="lambda.min")
hc[1:5]
```
    [1] -0.20607811 -0.11598244  0.00000000 -0.02551350 -0.06373397

## type.measure = "cindex"
```r
sfit = cv_rpair(x, S, loss_type="sqh", nlambda=25, type.measure="cindex", alignment="fraction")
plot(sfit, ggplot=F)
```

<img src="imgs/cv_sfit_plot.png" width="665" height="455" />

```r
spred = predict(sfit, newx = x, s="lambda.1se")
spred[1:5]
```
        [1]  1.7841666 -2.0529450 -0.7171475 -0.9923465 -0.3054223
    
```r
sc = coef(sfit, s="lambda.min")
sc[1:5]
```
    [1] -0.5366175 -0.4487022 -0.2064157 -0.3719257 -0.3550080

```r
lfit = cv_rpair(x, S, loss_type="log", nlambda=25, type.measure="cindex", alignment="fraction")
plot(lfit)
```

<img src="imgs/cv_lfit_plot.png" width="665" height="455" />

```r
lpred = predict(lfit, newx = x, s="lambda.min")
lpred[1:5]
```
    [1]  12.405622 -15.946027  -6.851173  -5.647310  -2.843435
    
```r
lc = coef(lfit, s="lambda.1se")
lc[1:5]
```
    [1] -1.1123256 -0.9166612 -0.2898089 -0.6603276 -0.6496257
