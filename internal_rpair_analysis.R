rm(list = ls())

# USER INPUT -----------------------------------
# model parameters
nreps = 100
nlambda = 100
g_type_measure = "deviance"  # type.measure for glmnet
r_type_measure = "deviance"  # type.measure for rpair
r_alignment = "fraction"  # alignment for rpair
s = "lambda.min"  # values of penalty parameter lambda

# script parameters
num_clust = 20
rand_seed = 42
mfrow = c(3,5)

# file names
workdir = "~/rpair_analysis/"
infile = "rppa_nods"
outfile = "internal_rpair_analysis.rds"


# FUNCTION DEFINITIONS -----------------------------------
run_parallel_analysis <- function(poc, cl){
  print("----------------------------------------------")
  cat("\n", poc$data$ctype[1],"- \n")
  
  # set up the data
  x = as.matrix( scale( poc$data$rppa ) )
  S = poc$data$S
  
  # CV folds for each tr
  spids <- poc$spids
  rm(poc)
  
  clusterExport(cl, c("x", "S", "nlambda", "g_type_measure", "r_type_measure", "r_alignment", "s"), envir = environment())
  parApply(cl, spids, 2, function(trinds, alpha = 1){
    #apply(spids, 2, function(trinds, alpha = 0){
    cat("m. ")
    tryCatch({
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
      
      # train the four models: cox (glmnet), exponential (rpair), logistic (rpair), and squared hinge (rpair)
      models<- list( 
        gcv = tryCatch({cv.glmnet(xtr,Str, family = "cox", type.measure = g_type_measure, nlambda = nlambda, foldid = fids, alpha = alpha, pmax = pmx)}),
        ecv = tryCatch({cv_rpair(xtr, Str, loss_type = "exp", nlambda = nlambda, alignment = r_alignment, type.measure = r_type_measure, foldid = fids, alpha = alpha, pmax = pmx)}),
        lcv = tryCatch({cv_rpair(xtr, Str, loss_type = "log", nlambda = nlambda, alignment = r_alignment, type.measure = r_type_measure, foldid = fids, alpha = alpha, pmax = pmx)}),
        hcv = tryCatch(({cv_rpair(xtr, Str, loss_type = "sqh", nlambda = nlambda, alignment = r_alignment, type.measure = r_type_measure, foldid = fids, pmax = pmx)}))
      )
      rm(xtr, Str)
      
      list( 
        
        res = tryCatch({sapply( models, function(cv){
          # check if it is an error values returend are NA
          c(c_index.min = unname( survConcordance( Sts ~ predict(cv, xts, s = s) )$concordance ),
            c_index.1se = unname( survConcordance( Sts ~ predict(cv, xts) )$concordance ), 
            n_feature.min = sum(coef(cv,s = s) != 0 ),
            n_feature.1se = sum(coef(cv) != 0 ))})
        }),
        models = models
      )
    })
  })
  
  
}

fcvid <- function(S,k){
  if(sum(S[,2])<5) stop("not enough events!")
  
  # separate events and censored events
  inds1 = S[,2]==1
  inds0 = S[,2]==0
  
  # randomly assign events and censored events to each fold
  k = sample(seq(k))
  ids = rep(0, length(S[,2]))
  ids[inds1] = sample(rep(k, length = sum(inds1)))
  ids[inds0] = sample(rep(k, length = sum(inds0)))
  ids
}


# MAIN -----------------------------------

library(survival)
library(magrittr)
library(glmnet)
library(rpair)
library(parallel)

setwd(workdir)
logfile = paste("x", format(Sys.time(),"%Y%m%d"), "rppalog", sep = "_")

# check that outfile doesn't exist
if(file.exists(outfile)) stop(glue::glue("File {outfile} already exists."))

# initiate cluster
cl <- makeCluster(num_clust, outfile= logfile)
clusterExport(cl, as.vector( lsf.str() ))
clusterEvalQ(cl, library(survival))
clusterEvalQ(cl, library(glmnet))
clusterEvalQ(cl, library(rpair))

load(infile)
# remove missing values
dfnods = dfnods %>% lapply( na.omit)
# train test splits and inner cvs
dfnods<- dfnods %>% lapply( function(nod, nboo=nreps) 
  list( data = nod, spids = replicate(nboo, {
    ids = fcvid(nod$S, 5) < 5
    ids[ids] = fcvid(nod$S[ids,], 5)
    ids
  })))

# run analysis
set.seed(rand_seed)
re <- lapply(dfnods, run_parallel_analysis, cl=cl)

# visualize results as boxplots
par(mfrow = mfrow)
lapply(names(re), function(xx) re[[xx]] %>% sapply( function(x) x$res[1,]) %>% t %>% boxplot(main=xx) )

# save results
saveRDS(re, file = outfile)
