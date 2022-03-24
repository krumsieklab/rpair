Introduction
============

The 'Regularized PAIrwise Ranking', or rpair, package is a convenient 
toolbox for conducting analyses on high-dimensional omics data. rpair combines
a general framework for casting complex outcomes into
pairwise learning problems with a penalty function allowing for sparse
solutions. This results in models that are not only able to fully
utilize such data, but are far more interpretable than traditional 
methods such as the Cox proportional hazard model. The package supports four
types of loss functions: logistic loss, exponential loss, squared hinge
loss and huberized hinge loss.

Reference
=========

Buyukozkan, Chetnik, and Krumsiek. \"Regularized PAIrwise Ranking survival analysis:
RPAIR\". (unpublished) 2022

Installation instructions
=========================

The rpair package can be installed using the following command:

```r
require(devtools)
devtools::install_github(repo=\"krumsieklab/rpair\", subdir=\"rpair\")
```

Getting Started
===============
```r
# generate survival data
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

library(rpair)

cv = cv_rpair(x, S, foldid=fids, nlambda=25, loss_type="exp", alpha = alpha,
                pmax = pmx, alignment = "fraction", keep = T)

plot(cv)

pr = predict(cv, x)
pr[1:5]
```

<img src="tutorials/imgs/get_start_plot.png" width="665" height="455" />

    [1]  0.79511265 -1.52421967 -0.53790436 -0.74291613 -0.27140260  0.01558794 -0.57333305  1.50476447  1.01263183
    [10] -1.13966754

Tutorials
=========

Detailed examples illustrating the full functionalities of the package
are provided in the following tutorials:

-   [Tutorial 1: the rpair function](https://github.com/krumsieklab/rpair/blob/master/tutorials/01_the_rpair_function.md)

-   [Tutorial 2: the cv_rpair function](https://github.com/krumsieklab/rpair/blob/master/tutorials/02_the_cv_rpair_function.md)

-   [Tutorial 3: analysis of dataset x](https://github.com/krumsieklab/rpair/blob/master/tutorials/03_analysis_of_dataset_x.md)
