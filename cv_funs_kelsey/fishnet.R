

srpoisson <- function(x, S, ...){
  xh = get_riskset(S) %*% x
  #cat(". ")
  glmnet( -xh, rep(0, nrow(xh)), family = "poisson", intercept = F, ...) 
} 

cv.srpoisson <- function(x, S, finds, 
                         pmx = NULL,
                         k = 33, # # of lambda to search
                         lfactor = ifelse(nobs < nvars, 0.01, 1e-04) ,
                         type.measure = "mse",
                         ... ){
  nvars = ncol(x)
  nobs = nrow(x)
  
  # get base object
  poi.obj = srpoisson(x, S, lambda.min.ratio = lfactor, pmax = pmx, nlambda = k, ...)
  lambda = poi.obj$lambda
  
  # get cv folds 
  fits <- lapply( seq(max(finds)), function(i) 
    srpoisson(x[finds != i,], S[finds != i,], lambda = lambda, ...) )
  
    
  #function(type.measure){
  cvstuff = mcv.fishnet(fits, poi.obj$lambda, x, S, foldid = finds, type.measure = type.measure)
  
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
  
  # I'm missing getmin()
  #lamin = getmin(lambda, cvm, cvsd)
  
  #obj = c(out, as.list(lamin))
  obj = out
  class(obj) = c( "cv.glmnet", "cv.srpoisson")
  obj
  #}
}

# calculate accuracy measure
mcv.fishnet <-
  function (outlist, lambda, x, y, weights=NULL, offset = NULL, foldid, type.measure = "deviance", 
            grouped = FALSE, keep = FALSE) 
  {
    
    #---
    
    cvs <-
      lapply(seq(max(foldid)), function(i){
        inds = foldid == i
        xh = get_riskset(y[inds,]) %*% x[inds,]
        ii = rep(i,nrow(xh))
        # print(nrow(xh))
        list(x=xh, i = ii)
      })
    x = do.call(rbind, lapply(cvs, `[[`,"x")) * (-1)
    y = x[,1]*0
    foldid = do.call(c, lapply(cvs, `[[`,"i"))
    rm(cvs)
    
    weights = y + 1
    #---
    
    if (!is.null(offset)) {
      is.offset = TRUE
      offset = drop(offset)
    }
    else is.offset = FALSE
    mlami = max(sapply(outlist, function(obj) min(obj$lambda)))
    which_lam = lambda >= mlami
    devi = function(y, eta) {
      deveta = y * eta - exp(eta)
      devy = y * log(y) - y
      devy[y == 0] = 0
      2 * (devy - deveta)
    }
    predmat = matrix(NA, length(y), length(lambda))
    nfolds = max(foldid)
    nlams = double(nfolds)
    for (i in seq(nfolds)) {
      which = foldid == i
      fitobj = outlist[[i]]
      if (is.offset) 
        off_sub = offset[which]
      preds = predict(fitobj, x[which, , drop = FALSE], s = lambda[which_lam], 
                      newoffset = off_sub)
      nlami = sum(which_lam)
      predmat[which, seq(nlami)] = preds
      nlams[i] = nlami
    }
    N = length(y) - apply(is.na(predmat), 2, sum)
    cvraw = switch(type.measure, mse = (y - exp(predmat))^2, 
                   mae = abs(y - exp(predmat)), deviance = devi(y, predmat))
    if ((length(y)/nfolds < 3) && grouped) {
      warning("Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold", 
              call. = FALSE)
      grouped = FALSE
    }
    if (grouped) {
      cvob = cvcompute(cvraw, weights, foldid, nlams)
      cvraw = cvob$cvraw
      weights = cvob$weights
      N = cvob$N
    }
    cvm = apply(cvraw, 2, weighted.mean, w = weights, na.rm = TRUE)
    cvsd = sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, weighted.mean, 
                      w = weights, na.rm = TRUE)/(N - 1))
    out = list(cvm = cvm, cvsd = cvsd, type.measure = type.measure)
    if (keep) 
      out$fit.preval = predmat
    out
  }




