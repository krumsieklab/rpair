
A_ <- function(X){
  # --- old version
  # y = rep(c(1,-1),nrow(X))[1:nrow(X)]
  # X = sweep(X,1,y,`*`)
  
  # new version
  X = rbind(X[1:2,],X)
  y = c(1,1,rep(0, nrow(X)-2))
  # ----------
  
  list(x = X, y = y)
}

srbinomial <- function(x, S, ...){
  aa = A_( get_riskset(S) %*% x )
  glmnet( aa$x, aa$y, family = "binomial", intercept = F, ...) 
} 

cv.srbinomial <- function(x, S, finds, 
                          pmx = 90,
                          k = 33, # # of lambda to search
                          lfactor = ifelse(nobs < nvars, 0.01, 1e-04) ,
                          type.measure = "deviance",
                          ... ){
  nvars = ncol(x)
  nobs = nrow(x)
  
  # get base object
  poi.obj = srbinomial(x, S, lambda.min.ratio = lfactor, pmax = pmx, nlambda = k, ...)
  lambda = poi.obj$lambda
  
  # get cv folds 
  fits <- lapply( seq(max(finds)), function(i) 
    srbinomial(x[finds != i,], S[finds != i,], lambda = lambda, ...) )
  
  
  cvstuff = mcv.lognet(fits, poi.obj$lambda, x, S, foldid = finds, type.measure = type.measure)
  
  nz = sapply(predict(poi.obj, type = "nonzero"), length)
  cvm = cvstuff$cvm
  cvsd = cvstuff$cvsd
  nas = is.na(cvsd)
  if (any(nas)) {
    lambda = lambda[!nas]
    cvm = cvm[!nas]
    cvsd = cvsd[!nas]
    nz = nz[!nas]
  }
  cvname = cvstuff$type.measure
  names(cvname) = cvstuff$type.measure
  
  out = list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvup = cvm + 
               cvsd, cvlo = cvm - cvsd, nzero = nz, name = cvname, glmnet.fit = poi.obj)
  
  lamin = getmin(lambda, cvm, cvsd)
  
  obj = c(out, as.list(lamin))
  class(obj) = c( "cv.glmnet", "cv.srbinomial")
  obj
}


mcv.lognet <-
  function (outlist, lambda, x, y, weights=NULL, offset = NULL, foldid, type.measure = "deviance", 
            grouped = FALSE, keep = FALSE) 
  {
    #---
    
    cvs <-
      lapply(seq(max(foldid)), function(i){
        inds = foldid == i
        aa = A_(get_riskset(y[inds,]) %*% x[inds,])
        ii = rep(i,nrow(aa$x))
        # print(nrow(xh))
        list(x=aa$x, y = aa$y ,i = ii)
      })
    x = do.call(rbind, lapply(cvs, `[[`,"x")) * (-1)
    y = do.call(c, lapply(cvs, `[[`,"y"))
    foldid = do.call(c, lapply(cvs, `[[`,"i"))
    rm(cvs)
    
    weights = y*0 + 1
    #---
    
    prob_min = 1e-05
    prob_max = 1 - prob_min
    nc = dim(y)
    if (is.null(nc)) {
      y = as.factor(y)
      ntab = table(y)
      nc = as.integer(length(ntab))
      y = diag(nc)[as.numeric(y), ]
    }
    N = nrow(y)
    nfolds = max(foldid)
    if ((N/nfolds < 10) && type.measure == "auc") {
      warning("Too few (< 10) observations per fold for type.measure='auc' in cv.lognet; changed to type.measure='deviance'. Alternatively, use smaller value for nfolds", 
              call. = FALSE)
      type.measure = cvtype("deviance", "lognet")
    }
    if ((N/nfolds < 3) && grouped) {
      warning("Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold", 
              call. = FALSE)
      grouped = FALSE
    }
    if (!is.null(offset)) {
      is.offset = TRUE
      offset = drop(offset)
    }
    else is.offset = FALSE
    mlami = max(sapply(outlist, function(obj) min(obj$lambda)))
    which_lam = lambda >= mlami
    predmat = matrix(NA, nrow(y), length(lambda))
    nlams = double(nfolds)
    for (i in seq(nfolds)) {
      which = foldid == i
      fitobj = outlist[[i]]
      if (is.offset) 
        off_sub = offset[which]
      preds = predict(fitobj, x[which, , drop = FALSE], s = lambda[which_lam], 
                      newoffset = off_sub, type = "response")
      nlami = sum(which_lam)
      predmat[which, seq(nlami)] = preds
      nlams[i] = nlami
    }
    if (type.measure == "auc") {
      cvraw = matrix(NA, nfolds, length(lambda))
      good = matrix(0, nfolds, length(lambda))
      for (i in seq(nfolds)) {
        good[i, seq(nlams[i])] = 1
        which = foldid == i
        for (j in seq(nlams[i])) {
          cvraw[i, j] = glmnet:::auc.mat(y[which, ], predmat[which, 
                                                    j], weights[which])
        }
      }
      N = apply(good, 2, sum)
      weights = tapply(weights, foldid, sum)
    }
    else {
      ywt = apply(y, 1, sum)
      y = y/ywt
      weights = weights * ywt
      N = nrow(y) - apply(is.na(predmat), 2, sum)
      cvraw = switch(type.measure, mse = (y[, 1] - (1 - predmat))^2 + 
                       (y[, 2] - predmat)^2, mae = abs(y[, 1] - (1 - predmat)) + 
                       abs(y[, 2] - predmat), deviance = {
                         predmat = pmin(pmax(predmat, prob_min), prob_max)
                         lp = y[, 1] * log(1 - predmat) + y[, 2] * log(predmat)
                         ly = log(y)
                         ly[y == 0] = 0
                         ly = drop((y * ly) %*% c(1, 1))
                         2 * (ly - lp)
                       }, class = y[, 1] * (predmat > 0.5) + y[, 2] * (predmat <= 
                                                                         0.5))
      if (grouped) {
        cvob = cvcompute(cvraw, weights, foldid, nlams)
        cvraw = cvob$cvraw
        weights = cvob$weights
        N = cvob$N
      }
    }
    cvm = apply(cvraw, 2, weighted.mean, w = weights, na.rm = TRUE)
    cvsd = sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, weighted.mean, 
                      w = weights, na.rm = TRUE)/(N - 1))
    out = list(cvm = cvm, cvsd = cvsd, type.measure = type.measure)
    if (keep) 
      out$fit.preval = predmat
    out
  }


