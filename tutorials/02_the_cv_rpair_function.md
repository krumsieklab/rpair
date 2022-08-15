Tutorial 2: the cv_rpair function
================

This tutorial demonstrates how to use the cv_rpair function and its utility functions. The cv_rpair function performs
k-fold cross-validation of the rpair function. The cv_rpair functions supports two kinds of input: survival and continuous. Like rpair, cv_rpair has three corresponding utility functions: predict, coef, and plot. 
Two important parameters featured in this tutorial are type.measure. The type.measure parameter allows  the user to select which loss function to use for cross-validation - "deviance" (default) or "cindex". Selecting "deviance" uses the selected loss function itself (e.g. exponential loss) while selecting "cindex" uses the concordance measure.

``` r
library(rpair)
library(magrittr)
library(survival)
```

## Random Survival Dataset and Fold IDs
```r

# data and outcome
x = rpair::ds2_x
y = rpair::ds2_y

# fold ids
set.seed(42)
fids <- sample(4, 200, T)

```

```r
# first 4 rows and columns of the data matrix
x[1:4,1:4]
```
               [,1]       [,2]       [,3]        [,4]
    [1,] -0.7943683 -0.7508653 -0.4146978  0.02697822
    [2,]  0.1972575 -0.5664613 -0.3503992  2.62760469
    [3,]  1.0017043  0.6997741  0.3240541 -0.99577921
    [4,]  1.2888254 -0.6372901  0.1341382  0.39321401


```r
# firs five rows and columns of the survival outcome matrix
y[1:5,]
```
         time status
    [1,]   54      1
    [2,]   59      1
    [3,]   52      1
    [4,]  141      0
    [5,]  200      1

## type.measure = "deviance"
```r
efit = cv_rpair(x, y, foldid=fids, nlambda=25, type_measure = "deviance", alignment = "fraction", 
                keep = T, loss_type="exp")
plot(efit)
```

<img src="imgs/cv_efit_plot.png" width="665" height="455" />

```r
epred = predict(efit, newx = x, s="lambda.min")
epred[1:5]
```
    [1]  0.36584855 -0.03888895 -0.06310045 -0.46314332 -1.02546453
    
```r
ec = coef(efit, s="lambda.1se")
ec[104:108]
```
    [1] -0.03163990  0.00000000 -0.01810281  0.00000000  0.03404023

```r
hfit = cv_rpair(x, y, loss_type="huh", foldid=fids, nlambda=25, type_measure = "deviance",
                keep = T, alignment="fraction", delta=3)
plot(hfit, ggplot=F)
```

<img src="imgs/cv_hfit_plot.png" width="665" height="455" />

```r
hpred = predict(hfit, x, s="lambda.1se")
hpred[1:5]
```
    [1]  0.26900400 -0.03424010 -0.05550277 -0.31228996 -0.70522639

```r
hc = coef(hfit, s="lambda.min")
hc[1:5]
```
    [1] 0.04802845 0.00000000 0.00000000 0.03425251 0.00000000

## type.measure = "cindex"
```r
sfit = cv_rpair(x, y, loss_type="sqh", nlambda=25, type_measure="cindex", alignment="fraction")
plot(sfit, ggplot=F)
```

<img src="imgs/cv_sfit_plot.png" width="665" height="455" />

```r
spred = predict(sfit, newx = x, s="lambda.min")
spred[1:5]
```
        [1]  0.12931793 -0.09054778 -0.07009912 -0.09598913 -0.41587683
    
```r
sc = coef(sfit, s="lambda.min")
sc[20:24]
```
    [1]  0.000000000  0.000000000 -0.012451232 -0.009717556  0.000000000

```r
lfit = cv_rpair(x, y, loss_type="log", nlambda=25, type_measure="cindex", alignment="fraction")
plot(lfit)
```

<img src="imgs/cv_lfit_plot.png" width="665" height="455" />

```r
lpred = predict(lfit, newx = x, s="lambda.min")
lpred[1:5]
```
    [1]  0.12765055 -0.12895854 -0.09370281 -0.05837271 -0.60130044
    
```r
lc = coef(lfit, s="lambda.min")
lc[86:90]
```
    [[1]  0.00000000  0.00000000  0.00000000 -0.05780114  0.07670461
