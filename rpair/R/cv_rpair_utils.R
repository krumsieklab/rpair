# cross-validation random folds
get_folds <- function(y, nfolds){
  y <- cbind(index = 1:length(y), y)
  sorted <- y[order(y[,2]),]
  fold_df <- create_fold_splits(sorted, nfolds)

  # handle indices with unassigned folds
  nas <- sum(is.na(fold_df$fold))
  if(nas <= nfolds){
    fold_df[is.na(fold_df$fold),"fold"] <- sample(seq(1:nfolds))[1:nas]
  }else{
    fold_df[is.na(fold_df$fold),"fold"] <- sample(seq(1:nfolds), size=nas, replace = T)
  }
  folds <- merge(sorted, fold_df, by="index")
  folds
}

# cross-validation stratified folds
get_stratified_folds <- function(S, nfolds){
  Smat <- as.data.frame(as.matrix(S))
  Smat <- cbind(index = 1:dim(Smat)[1], Smat)

  # sort survial outcome by time variable
  sorted <- Smat[order(Smat[,2]),]
  eventlist <- split(sorted, f=sorted[,3])
  fold_df_list <- lapply(eventlist, create_fold_splits, nfolds)
  fold_df <- do.call(rbind.data.frame,fold_df_list)

  # handle indices with unassigned folds
  nas <- sum(is.na(fold_df$fold))
  if(nas <= nfolds){
    fold_df[is.na(fold_df$fold),"fold"] <- sample(seq(1:nfolds))[1:nas]
  }else{
    fold_df[is.na(fold_df$fold),"fold"] <- sample(seq(1:nfolds), size=nas, replace = T)
  }
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
                                 fishnet = "Exponential Loss",
                                 phuhnet = "Huberized Hinge Loss",
                                 psqhnet = "Squared Hinge Loss")
  }else{
    names(type.measure) <- "C-Index"
  }
  type.measure
}

getopt_cv_rpair <- function (lambda, cvm, cvsd, cvname)
{
  cvm = cvm[!is.na(lambda)]
  cvsd = cvsd[!is.na(lambda)]
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
huh_loss <- function(x, delta) ((1 - -1*x - 0.5 * delta) * (-1*x <= 1 - delta) + 0.5 * (1 - -1*x)^2/delta * (-1*x <= 1) * (-1*x > 1 - delta) )
sqh_loss <- function(x) ((x >= 0) * (1+x)^2 )

cv_deviance <- function(predmat, y, loss_type, fids, delta, use_houw){
  ffs <- list(exp = exp_loss, log = log_loss, huh = huh_loss, sqh = sqh_loss)
  pi_all = y_to_pairs(y)

  lo =
    sapply(seq(unique(fids)), function(i){
      if(use_houw){
        # training set samples
        j = which(fids != i)
        # all_pairs - training_set_pairs
        pir =  pi_all[!( (pi_all[,1] %in% j) & (pi_all[,2] %in% j) ),]
      }else{
        # test set samples
        j = which(fids == i)
        # test_pairs only
        pir = pi_all[( (pi_all[,1] %in% j) & (pi_all[,2] %in% j) ),]
      }

      # loss margin z i.e. loss_function(z)
      Z = predmat[[i]][pir[,1],] -  predmat[[i]][pir[,2],]
      # Z is a matrix of size:n_valid_pairs X n_lambda

      # pick loss function based on the model loss type
      ff = ffs[[loss_type]] # ffs is the list of functions

      # for each lambda calculate deviance for fold i
      if(loss_type != "phuhnet"){
        colMeans( ff(Z), na.rm = T )
      }else{
        colMeans( ff(Z, delta), na.rm = T )
      }
    })

  lo

}

cv_concordance <- function(predmat, y, foldid, use_houw){
  Yh = do.call(predmat, what = rbind)
  # folds, not fold ids
  z = unlist( lapply(seq(predmat), function(i) rep(i, nrow(predmat[[i]]))) )
  # inds of training sets
  ii = unlist(lapply(seq(predmat), function(i) foldid != i))
  # S for each fold
  Sf = rep(y, length(predmat))
  if(use_houw){
    cvfl = apply(Yh, 2, function(yh)
      tryCatch(calc_conc(Sf, yh, ii, z)["houw",], error = function(er) c(NA,NA)))
  }else{
    cvfl = apply(Yh, 2, function(yh)
      tryCatch(calc_conc(Sf, yh, ii, z)["test",], error = function(er) c(NA,NA)))
  }
  cvfl
}

calc_conc <- function(S, yh, ii, folds){

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
    d = c(full = 1-d0, train = 1-d1, houw = 1-dh, test = 1-d2$concordance),
    v = c(full = v0, train = v1, houw = vh, test = d2$var )
  )
}

cv_stats <- function (cvfl, lambda, nz, conc=F)
{
  if(conc){
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
