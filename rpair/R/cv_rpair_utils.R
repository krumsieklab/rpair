# cross-validation stratified folds
get_stratified_folds <- function(S, nfolds){
  Smat <- as.data.frame(as.matrix(S))
  Smat <- cbind(index = 1:dim(Smat)[1], Smat)
  sorted <- Smat[order(Smat$time),]
  eventlist <- split(sorted, f=sorted$status)
  fold_df_list <- lapply(eventlist, create_fold_splits, nfolds)
  fold_df <- do.call(rbind.data.frame,fold_df_list)
  fold_df[is.na(fold_df$fold),"fold"] <- sample(seq(1:nfolds))[1:sum(is.na(fold_df$fold))]
  strat_folds <- merge(sorted, fold_df, by="index")
  strat_folds
}

create_fold_splits <- function(df, nfolds){
  d <- df[,"index"]
  bins <- split(d, ceiling(seq_along(d)/nfolds))
  nbins = length(bins)
  fold_list <- lapply(bins, function(x){
    if(length(x)==nfolds){
      cbind(index=x, fold=sample(seq(1,nfolds)))
    }else{
      cbind(index=x, fold=rep(NA,length(x)))
    }
  })
  folds <- do.call(rbind,fold_list)
  folds
}



rpair_cvtype <- function(type.measure, losstype){
  if(type.measure=="deviance"){
    type.measure = losstype
    names(type.measure) = switch(losstype, lognet = "Log Loss",
                                 fishnet = "Exponential Loss")
  }else{
    names(type.measure) <- "C-Index"
  }
  type.measure
}

getopt_cv_rpair <- function (lambda, cvm, cvsd, cvname)
{
  if (match(cvname, "cindex", 0)) cvm = -cvm
  cvmin = min(cvm, na.rm = TRUE)
  idmin = cvm <= cvmin
  lambda.min = max(lambda[idmin], na.rm = TRUE)
  idmin = match(lambda.min, lambda)
  semin = (cvm + cvsd)[idmin]
  id1se = cvm <= semin
  lambda.1se = max(lambda[id1se], na.rm = TRUE)
  id1se = match(lambda.1se, lambda)
  index = matrix(c(idmin, id1se), 2, 1, dimnames = list(c("min",
                                                          "1se"), "Lambda"))
  list(lambda.min = lambda.min, lambda.1se = lambda.1se, index = index)
}


rpair_buildPredmat <- function(outlist, nlambda, lambda, x, foldid, alignment){
  if(alignment=="lambda") nlambda = length(lambda)
  pred_list = list()
  nfolds = max(foldid)
  nlams = double(nfolds)
  for (i in seq(nfolds)) {
    which = foldid == i
    fitobj = outlist[[i]]
    predmat = matrix(NA,nrow(x),nlambda)
    preds = switch(alignment, fraction = predict(fitobj, x),
                   lambda = predict(fitobj, x, s = lambda))
    nlami = min(ncol(preds), nlambda)
    predmat[, seq(nlami)] = preds[, seq(nlami)]
    rn = rownames(x)
    sn = paste("s", seq(0, length = nlambda), sep = "")
    dimnames(predmat) = list(rn, sn)
    pred_list[[i]] = predmat
  }
  pred_list
}

exp_loss <- function(x) exp(x)
log_loss <- function(x) log(1+exp(x))
#conc_loss <- function(x) 1*(x<0)

cv_loss <- function(predmat, pi_all, type_measure, fids){
  ffs <- list(fishnet = exp_loss, lognet = log_loss)#, cindex = conc_loss)

  lo =
    sapply(seq(unique(fids)), function(i){
      # fold i
      # training set samples
      j = which(fids != i)

      # all_pairs - training_set_pairs
      pir =  pi_all[!( (pi_all[,1] %in% j) & (pi_all[,2] %in% j) ),]
      npairs = nrow(pir)
      # KC: do we want to keep this number?

      # loss margin z i.e. loss_function(z)
      Z = predmat[[i]][pir[,1],] -  predmat[[i]][pir[,2],]
      # Z is a matrix of size:n_valid_pairs X n_lambda

      # pick loss function based on the model loss type
      ff = ffs[[type_measure]] # ffs is the list of functions

      # for each lambda calculate deviance for fold i
      # colMeans(exp_loss(Z),na.rm = T)
      colMeans( ff(Z), na.rm = T )
    })

  lo

}

cv_stats <- function (cvfl, lambda, nz, houw=F)
{
  if(houw){
    cvm = cvfl[1,]
    cvsd = sqrt(cvfl[2,])
  }else{
    cvm = rowMeans(cvfl, na.rm = T)
    cvsd = apply(cvfl, 1, sd, na.rm = T) / sqrt(ncol(cvfl))
  }
  nas = is.na(cvsd)
  if (any(nas)) {
    lambda = lambda[!nas]
    cvm = cvm[!nas]
    cvsd = cvsd[!nas]
    nz = nz[!nas]
  }
  list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvup = cvm +
         cvsd, cvlo = cvm - cvsd, nzero = nz)
}

cv_houw_loss <- function(S, yh, ii, folds){

  df = data.frame(S=S, yh=yh, z=folds)
  # concordance on full data
  d0 = concordance(S~yh+strata(z), df, keepstrata = T, reverse = T)
  d0$count = rbind(d0$count)
  # concordance on train data
  d1 = concordance(S~yh+strata(z), df[ii,], keepstrata = T, reverse = T)
  d1$count = rbind(d1$count)

  #estimated concordance for houw pairs
  m = colSums(d0$count) - colSums(d1$count)
  dh = unname( (m[1]+m[3]/2)/sum(m[1:3]) )

  #--- estimate se of houw pairs
  # number of all pairs and its variance
  n0 = sum(d0$count[,1:3])
  v0 = d0$var
  d0 = d0$concordance

  # number of train pairs and variance
  n1 = sum(d1$count[,1:3])
  v1 = d1$var
  d1 = d1$concordance

  # number of houw pairs
  nh =  sum(m[1:3])
  # variance of h pairs
  vh = unname( ( n0*v0*(n1+n0) -n1*v1*n1 - n1*(d0-d1)^2 - nh*(d0-dh) )/nh )/nh

  d2 =  concordance(S~yh+strata(z), df[!ii,], reverse = T)
  cbind(
    d = c(full = d0, train = d1, houw = dh, test = d2$concordance),
    v = c(full = v0, train = v1, houw = vh, test = d2$var )
  )
}
