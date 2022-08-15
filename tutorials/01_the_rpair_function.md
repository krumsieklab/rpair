Tutorial 1: the rpair function
================

This tutorial demonstrates how to use the rpair function and its utility functions.  The rpair function supports four types of loss functions: exponential (default), logistic, squared hinge, and huberized hinge loss. Each loss type is demonstrated below along with one or more of rpair's utility functions.

NOTE: The rpair function performs model fitting on a dataset, however, in most cases it is desirable to optimize the models parameters using cross-validation. The cv_rpair function performs k-fold cross-validation for the rpair function. In most use cases, the user will call cv_rpair rather than call the rpair function directly.  The tutorial for the cv_rpair function can be found [here](https://github.com/krumsieklab/rpair/blob/master/tutorials/02_the_cv_rpair_function.md).


``` r
library(rpair)
library(magrittr)
library(survival)
```

## Example Survival Dataset
```r
# load example data
x = rpair::ds2_x
y = rpair::ds2_y
```

```r
# first four rows and columns of the data matrix
x[1:4,1:4]
```
               [,1]       [,2]       [,3]        [,4]
    [1,] -0.7943683 -0.7508653 -0.4146978  0.02697822
    [2,]  0.1972575 -0.5664613 -0.3503992  2.62760469
    [3,]  1.0017043  0.6997741  0.3240541 -0.99577921
    [4,]  1.2888254 -0.6372901  0.1341382  0.39321401 

```r
# first five rows of survival data
y[1:5,1:2]
```
         time status
    [1,]   54      1
    [2,]   59      1
    [3,]   52      1
    [4,]  141      0
    [5,]  200      1



## Fit with exponential Loss
```r
efit = rpair(x, y, loss_type="exp", nlambda=25, pmax=50)
# generate trace plot
plot(efit)
```

<img src="imgs/efit_plot.png" width="665" height="455" />

```r
# extract coefficients from rpair object
ec <- coef(efit)
ec[22:24,10:14]
```
                 s9         s10         s11          s12           s13
    V22 -0.02363581 -0.01516192 -0.00585423 -0.001641664 -0.0006924081
    V23 -0.11643732 -0.12279514 -0.12787731 -0.131600482 -0.1317539180
    V24 -0.04682293 -0.05965208 -0.07098242 -0.079431407 -0.0845784596


## Fit with logistic Loss
```r
lfit = rpair(x, y, loss_type="log", nlambda=25, pmax=50)
# generate trace plot
plot(lfit, xvar="dev")
```

<img src="imgs/lfit_plot.png" width="665" height="455" />

```r
# predict values from rpair object
lp <- predict(lfit, newx=x)
lp[1:4, 2:5]  # skip the intercept column
```
                  s1          s2         s3         s4
    [1,]  0.03509707  0.12765055  0.2677591  0.4346283
    [2,] -0.05504652 -0.12895854 -0.1974425 -0.1811280
    [3,] -0.05865191 -0.09370281 -0.1365903 -0.1439619
    [4,] -0.02853772 -0.05837271 -0.1900639 -0.3370274
    
    
## Fit with squared Hinge Loss
```r
sfit = rpair(x, y, loss_type="sqh", nlambda=25, pmax=50)
plot(sfit, xvar="lambda")
```

<img src="imgs/sfit_plot.png" width="665" height="455" />

```r
sn <- predict(sfit, type = "nonzero")
sn$s2 # [JK it seems like this fit is inspected different than the other two above and the fourth one below... why is that?]
```
    [1]   9  22  84  89  90  95 104 108 118 142 143 173 175 198
```r
sn$s8
```
    [1]   1   4   9  13  14  15  16  18  20  22  23  24  26  27  35  37  38  45  46  47  48  49  51  52  55  57  60  74     75  78
    [31]  80  84  86  88  89  90  91  94  95  96  97  99 100 101 102 103 104 106 108 112 115 117 118 121 122 124 125 126     129 136
    [61] 137 139 140 142 143 145 148 154 161 163 165 167 173 175 176 180 181 184 188 190 191 197 198  


## Huberized Hinge Loss
```r
hfit = rpair(x, y, loss_type="huh", nlambda=25, pmax=50)
plot(hfit, xvar="norm")
```

<img src="imgs/hfit_plot.png" width="665" height="455" />

```r
hc <- predict(hfit, type="coefficient")
hc[21:23,10:14]
```

                 s9          s10          s11          s12         s13
    V21  0.00000000  0.000000000  0.000000000  0.003733044  0.01349161
    V22 -0.01286632 -0.008549847 -0.002994905  0.000000000  0.00000000
    V23 -0.06644389 -0.069723974 -0.070907255 -0.071466062 -0.07189343
