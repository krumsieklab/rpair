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
[JK: intro needs to be adapted, up to this point it could be glmnet. currently discussed in slack] 

Reference
=========

Buyukozkan, Chetnik, and Krumsiek. \"Regularized PAIrwise Ranking survival analysis:
RPAIR\". (unpublished) 2022

Installation instructions
=========================

The rpair package can be installed using the following command:

```r
require(devtools)
devtools::install_github(repo="krumsieklab/rpair", subdir="rpair")
```

Getting Started
===============
```r

# libraries
library(rpair)
library(magrittr)
library(survival)

set.seed(42)

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

fids = sample(4, ks, T) # [JK, can this be deleted?]

# [JK, @mubu, I think you said this can be deleted?]
alpha = 1
# the maximum number of variables ever to be nonzero
pmx = min(ncol(X), sum(S[,2]))
if(alpha < 0.5 ) pmx = ncol(X)+1


# run cross-validated rpair
cv = cv_rpair(X, S, foldid=fids, nlambda=25, loss_type="exp", alpha = alpha,
                pmax = pmx, alignment = "fraction", keep = T)

# lambda plot
plot(cv)

# apply model to input data
pr = predict(cv, X)
# show predicted survival times
pr[1:5]
```

[JK, @Kelsey, it looks like the vector below is not actually an output of the code above. Code is 1:5, but below is 10 entries]

<img src="tutorials/imgs/get_start_plot.png" width="665" height="455" />

    [1]  0.79511265 -1.52421967 -0.53790436 -0.74291613 -0.27140260  0.01558794 -0.57333305  1.50476447  1.01263183
    [10] -1.13966754



Tutorials
=========

Detailed examples illustrating the full functionality of the package
are provided in the following tutorials:

-   [Tutorial 1: the rpair function](https://github.com/krumsieklab/rpair/blob/master/tutorials/01_the_rpair_function.md)

-   [Tutorial 2: the cv_rpair function](https://github.com/krumsieklab/rpair/blob/master/tutorials/02_the_cv_rpair_function.md)

-   [Tutorial 3: analysis of dataset x](https://github.com/krumsieklab/rpair/blob/master/tutorials/03_analysis_of_dataset_x.md)
