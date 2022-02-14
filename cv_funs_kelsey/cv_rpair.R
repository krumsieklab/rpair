
# Change the names from '.' to '_'
# also support the case where no fold-ids are given
#   after finishing up
#   cindex - call concordance fun (already implemented)
#   loss - do in R; later consider moving to C
# add support for pairs passed to this function


#' Cross-Validation function for rpair_gloss
#' 
#' Performs k-fold cross-validation for rpair_gloss.
#' 
#' @param x Input matrix as in rpair_gloss.
#' @param y Response as in rpair_gloss.
#' @param nlambda Number of lambda values. Default: 100.
#' @param type.measure Loss to use for cross-validation. Available loss functions are "deviance" and "cindex".
#' @param nfolds Number of folds. Default: 10.
#' @param alignment Determines 
#' @param grouped Not implemented in this version. Default: FALSE.
#' @param keep If keep=TRUE, returns a list of fitted values for each fold and a corresponding list of foldids.
#'    Default: FALSE.
#' @param \dots Additional parameters to pass to the rpair_gloss function.
#' 
#' @return An object of class "cv_rpair" containing the following list of values:
#' 


cv_rpair <- function(x,
                     y,
                     lambda=NULL,
                     nlambda=100,
                     #type.measure = c("default", "mse", "deviance", "auc", "mae"),
                     #type.measure = c("deviance", "cindex"),
                     type.measure = c("deviance"),
                     nfolds = 10,
                     foldid = NULL,
                     alignment = c("lambda", "fraction"),
                     grouped = FALSE,
                     keep = FALSE,
                     ...
){

  
  type.measure = match.arg(type.measure)
  alignment = match.arg(alignment)
  if (!is.null(lambda) && length(lambda) < 2) 
    stop("Need more than one value of lambda for cv_rpair")
  if (!is.null(lambda) && alignment == "fraction") {
    warning("fraction of path alignment not available if lambda given as argument; switched to alignment=`lambda`")
    alignment = "lambda"
  }
  N = nrow(x)
  y = drop(y)
  cv.call = rpair.call = match.call(expand.dots = TRUE)
  which = match(c("type.measure", "nfolds", "foldid", "grouped", 
                  "keep"), names(rpair.call), FALSE)
  if (any(which)) 
    rpair.call = rpair.call[-which]
  rpair.call[[1]] = as.name("rpair")
  
  if (is.null(foldid)) 
    foldid = sample(rep(seq(nfolds), length = N))
  else nfolds = max(foldid)
  if (nfolds < 3) 
    stop("nfolds must be bigger than 3; nfolds=10 recommended")
  cv_rpair_raw(x, y, lambda, nlambda, type.measure, 
                   nfolds, foldid, alignment, grouped, keep, 
                   rpair.call, cv.call, ...)
  
}

cv_rpair_raw <- function(x, y, lambda, nlambda, type.measure, nfolds, 
                         foldid, alignment, grouped, keep, rpair.call, 
                         cv.call, ...){

  rpair.object = rpair_gloss(x, y, lambda = lambda, nlambda = nlambda, ...)
  rpair.object$call = rpair.call
  losstype = class(rpair.object)[[1]]
  #type.measure = rpair.cvtype(type.measure, losstype)
  
  N = nrow(x)
  
  if(alignment=="lambda") lambda = rpair.object$lambda
  nz = sapply(predict(rpair.object, type = "nonzero"), 
              length)
  outlist <- lapply( seq(max(foldid)), function(i) 
    rpair_gloss(x[foldid != i,], y[foldid != i,], lambda = lambda, nlambda = nlambda, ...))
  
  
  lambda = rpair.object$lambda
  class(outlist) = paste0(losstype, "list")
  predmat = rpair_buildPredmat(outlist, nlambda, lambda, x, foldid, 
                         alignment, y = y, grouped = grouped, 
                         type.measure = type.measure, family = class(rpair.object)[1])

  pi_all = y_to_pairs(y)
  cvfl = list(N = nrow(pi_all))
  cvfl$lo = cv_loss(predmat, pi_all, losstype, foldid)
  #out = cv_stats(cvfl, foldid, nfolds, lambda, nz)
  out = cv_stats(cvfl, lambda, nz)
  
  out = c(out, list(call = cv.call, rpair.fit = rpair.object))
  if (keep) 
    out = c(out, list(fit.preval = predmat, foldid = foldid))

  out
  
}

# I edited your implementation by deleting all .[x,], 
# becasue I need prediction for all samples to be able to calculate Houwelingen cv
rpair_buildPredmat <- function(outlist, nlambda, lambda, x, foldid, alignment, ...){
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

cv_loss <- function(predmat, pi_all, loss_type, fids){
  ffs <- list(fishnet = exp_loss, lognet = log_loss)
  
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
      ff = ffs[[loss_type]] # ffs is the list of functions 
      
      # for each lambda calculate deviance for fold i 
      # colMeans(exp_loss(Z),na.rm = T)
      colMeans( ff(Z), na.rm = T )
    })
  
  t(lo)
  
}

cv_stats <- function (cvfl, lambda, nz) #foldid, nfolds, lambda, nz) # add foldid and nfolds back later?
{
  # KC: ignore for now
  # if (grouped) {
  #   nlams = rep(dim(cvstuff$cvraw)[2], nfolds)
  #   cvstuff = cvcompute(cvstuff, foldid, nlams)
  # }
  cvm = with(cvfl, apply(lo, 2, mean, na.rm = TRUE))
  cvsd = with(cvfl, sqrt(apply(scale(lo, cvm, FALSE)^2, 
                                  2, mean, na.rm = TRUE)/(N - 1)))
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
  

