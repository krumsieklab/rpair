Introduction
============

The 'Regularized PAIrwise Ranking survival analysis', or rpair, package
is a convenient toolbox for conducting survival analyses on
high-dimensional omics data with censored survival times. rpair combines
a general framework for casting the survival prediction problem to a
pairwise learning problem with a penalty function allowing for sparse
solutions. This results in models that are not only able to fully
utilize censored data, but are far more interpretable than the
traditional Cox proportional hazard model. The package supports four
types of loss functions: logistic loss, exponential loss, squared hinge
loss and huberized hinge loss.

Reference
=========

Buyukozkan, et al. \"Regularized PAIrwise Ranking survival analysis:
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

library(rpair)

cv = cv_rpair(data, ...)

plot(cv)
pr = predict(cv, new_data)

# show plot and predictions
```

Tutorials
=========

Detailed examples illustrating the full functionalities of the package
are provided in the following tutorials:

-   [Tutorial 1: the rpair function](https://github.com/krumsieklab/rpair/blob/master/tutorials/01_the_rpair_function.md)

-   [Tutorial 2: the cv_rpair function](https://github.com/krumsieklab/rpair/blob/master/tutorials/02_the_cv_rpair_function.md)

-   Tutorial 3: example analysis function (?)
