Tutorial 3: analysis of dataset GSE3493NA
================

This tutorial demonstrates how to use rpair for performing a simple external analysis. The example analysis is performed using two public transcriptomics breast cancer datasets from the mcsurvdata package: GSE2990NA (model fitting) and GSE3493NA (model evaluation).

``` r
library(rpair)
library(magrittr)
library(survival)
library(mcsurvdata)
```

## Download and Tidy Datasets
```r

eh <- ExperimentHub()
dat <- query(eh, "mcsurvdata")
nda.brca <- dat[["EH1497"]]

# survival outcome 
S =  Surv(nda.brca$tev+1, nda.brca$evn == 1)
nda.brca = nda.brca[,!is.na(S)] 

# gene expression data
mrna = Biobase::assayData(nda.brca)$exprs %>% t

# survival outcome 
S =  data.frame(t = nda.brca$tev+1, e = nda.brca$evn == 1)

# cohorts 
Z = paste0(nda.brca$dataset, nda.brca$cohort)
cohorts = c("GSE2990NA","GSE3494NA")

nods <-
  lapply( structure( cohorts, names = cohorts ), function(i){
    print(i)
    inds = (Z==i)
    x = scale(mrna[inds,])
    y = Surv( S[inds,"t"], S[inds,"e"] )
    list(x=x, y=y)
  })
  

```

## Train on GSE2990NA
```r
# train data and outcome
x_tr = nods$GSE2990NA$x
y_tr = nods$GSE2990NA$y

# train model
set.seed(42)
fit = cv_rpair(y = y_tr, x = x_tr, loss_type = "log")

```

## Test on GSE3493NA
```r
# test data and outcome
x_te = nods$GSE3494NA$x
y_te = nods$GSE3494NA$y

# calculate concordance
conc = concordance(y_te~predict(fit, x_te, s = "lambda.min" ), reverse = T)$concordance

```
