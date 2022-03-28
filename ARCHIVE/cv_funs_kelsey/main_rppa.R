rm(list = ls())
nice <- function(){system("renice -n 19 -u buyukozk")}
system("pkill -u buyukozk -f slave")

library(survival)
library(magrittr)
library(glmnet)

setwd("~/survrank/")
source("fishnet.R")
source("lognet.R")
source("ssvm_utils.R")

library(parallel)
logfile = paste("x", format(Sys.time(),"%Y%m%d"), "rppalog", sep = "_")

# Initiate cluster
nice()
cl <- makeCluster(40, outfile= logfile)
clusterExport(cl, as.vector( lsf.str() ))
clusterEvalQ(cl, library(survival))
clusterEvalQ(cl, library(glmnet))

benchorse <- function(poc, cl){
  print("----------------------------------------------")
  cat("\n", poc$data$ctype[1],"- \n")
  
  # set up the data
  x = as.matrix( scale( poc$data$rppa ) )
  S = poc$data$S
  
  # CV folds for each tr
  spids <- poc$spids
  rm(poc)
  
  clusterExport(cl, c("x", "S"), envir = environment())
  parApply(cl, spids, 2, function(trinds, alpha = 1){
    #apply(spids, 2, function(trinds){
    cat("m. ")
    try({
      # weird glmnet error for surv.time == 0
      S[,1] = S[,1] + 1
      
      # training data 1,2,3,4,5:train data
      xtr = x[trinds>0, ]
      Str = S[trinds>0, ]
      
      # test data 0:test data
      Sts = S[trinds==0, ]
      xts = x[trinds==0, ]
      
      # cv ids
      fids = trinds[trinds>0] 
      
      # the maximum number of variables ever to be nonzero 
      pmx = min(ncol(x), sum(S[,2]))
      if(alpha < 0.5 ) pmx = ncol(x)+1
      
      models<- list( 
        gcv = cv.glmnet(xtr,Str, family = "cox", foldid = fids, alpha = alpha, pmax = pmx, nlambda = 25 ),
        pcv = cv.srpoisson(xtr, Str, fids,k = 25, type.measure = "deviance", alpha = alpha, pmx = pmx),
        bcv = cv.srbinomial(xtr, Str, fids, k = 25, type.measure = "auc", alpha = alpha, pmx = pmx)
      )
      rm(xtr, Str)
      
      list( 
        res = sapply( models, function(cv) 
          c(c_index.min = unname( survConcordance( Sts ~ predict(cv, xts, s = "lambda.min") )$concordance ),
            c_index.1se = unname( survConcordance( Sts ~ predict(cv, xts) )$concordance ), 
            n_feature.min = sum(coef(cv,s = "lambda.min") != 0 ),
            n_feature.1se = sum(coef(cv) != 0 )) ),
        models = models
      )
    })
  })
  
  
}

set.seed(42)
load("rppa_nods")
# remove missing values
dfnods = dfnods %>% lapply( na.omit)

# train test splits and inner cvs
dfnods<- dfnods %>% lapply( function(nod, nboo=100) 
  list( data = nod, spids = replicate(nboo, {
    ids = fcvid(nod$S, 5) < 5
    ids[ids] = fcvid(nod$S[ids,], 5)
    ids
  })))

re <- lapply(dfnods, benchorse, cl=cl)

par(mfrow = c(3,5))
lapply(names(re), function(xx) re[[xx]] %>% sapply( function(x) x$res[1,]) %>% {.[.==0]=NA;.} %>% t %>% boxplot(main=xx) )

save(file = "results_rppa_l1", rel1)

