rm( list = ls() )

out_file = "external_rpair_analysis_rerun.rds"

zeromat <- function (nvars, nalam, vnames, stepnames) {
  ca <- rep(0, nalam)
  ia <- seq(nalam + 1)
  ja <- rep(1, nalam)
  dd <- c(nvars, nalam)
  new("dgCMatrix", Dim = dd, Dimnames = list(vnames, stepnames),
      x = as.vector(ca), p = as.integer(ia - 1), i = as.integer(ja - 1))
}

library(mcsurvdata)
eh <- ExperimentHub()
dat <- query(eh, "mcsurvdata")
nda.brca <- dat[["EH1497"]]

# survival outcome 
S =  Surv(nda.brca$tev+1, nda.brca$evn == 1)
nda.brca = nda.brca[,!is.na(S)] 

# nda.crc <- dat[["EH1498"]]
nda.brca$ER.status %>% table
nda.brca$PGR.status %>% table
nda.brca$HER2.status %>% table

nda.brca$HER2.status %>% table( paste0(nda.brca$dataset, nda.brca$cohort))
# gene expression data
mrna = Biobase::assayData(nda.brca)$exprs %>% t

# survival outcome 
S =  data.frame(t = nda.brca$tev+1, e = nda.brca$evn == 1)

# cohorts 
Z = paste0(nda.brca$dataset, nda.brca$cohort)

# collect data 
nods <-
  lapply( structure( unique(Z), names = unique(Z) ), function(i){
    print(i)
    inds = (Z==i)
    x = scale(mrna[inds,])
    y = Surv( S[inds,"t"], S[inds,"e"] )
    list(x=x, y=y)
  })

# calculate cox ph results 
library(glmnet)
library(rpair)
methods = c("cox", "log", "exp", "hinge")

sub_nods <- nods[c("metabric1", "metabric2", "metabric3")]

run_extern_analysis_2 <- function(m, sub_nods){
  print(m)
  
  sub_seq = seq(sub_nods)
  res = lapply(sub_seq, function(i){
    print(glue::glue("Test: {names(sub_nods)[i]}"))
    # select test set
    x_te = sub_nods[[i]]$x
    y_te =  sub_nods[[i]]$y
    if(is.Surv(y_te)==FALSE) stop("y_te is not of class Surv!")
    
    # select and combine training sets
    ss = sub_seq[sub_seq != i]
    x_tr = rbind(sub_nods[[ss[1]]]$x, sub_nods[[ss[2]]]$x)
    y_tr = rbind(sub_nods[[ss[1]]]$y, sub_nods[[ss[2]]]$y)
    # convert to Surv class - weird error otherwise
    y_tr = Surv(y_tr[,1], y_tr[,2])
    if(is.Surv(y_tr)==FALSE) stop("y_tr is not of class Surv!")
    
    # cv glmnet 
    fit = switch(m, 
                 cox = cv.glmnet(y = y_tr, x = x_tr, nfolds = 5, family = "cox"),
                 log = cv_rpair(y = y_tr, x = x_tr, nfolds = 5, loss_type = "log", alignment = "fraction"),
                 exp = cv_rpair(y = y_tr, x = x_tr, nfolds = 5, loss_type = "exp", alignment = "fraction"),
                 hinge = cv_rpair(y = y_tr, x = x_tr, nfolds = 5, loss_type = "sqh", alignment = "fraction"))
    
    conc = concordance(y_te~predict(fit, x_te, s = "lambda.min" ), reverse = T)$concordance
    
    res = list("conc"=conc, "fit"=fit)
    res
    
  })
  
  res
}

start_time <- Sys.time()
res_list <- lapply(methods[1:3], run_extern_analysis_2, sub_nods)
end_time <- Sys.time()
end_time - start_time

saveRDS(res_list, file = out_file)
