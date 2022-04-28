#' Cross-Validation function for rpair_gloss and rpair_hinge.
#'
#' Performs k-fold cross-validation for rpair_gloss or rpair_hinge.
#'
#' @param x Input matrix as in rpair_gloss and rpair_hinge.
#' @param y Response as in rpair_gloss and rpair_hinge.
#' @param loss_type Loss function to use. One of c("exp", "log" ,"sqh", "huh"). Default: "exp".
#' @param nlambda Number of lambda values. Default: 100.
#' @param type.measure Loss to use for cross-validation. Available loss functions are "deviance" and "cindex".
#' @param nfolds Number of folds. Default: 5.
#' @param foldid Optional vector of values between 1 and nfolds identifying which fold each observation should be in.
#'   If provided, nofolds can be missing. Pairs caanot be provided as input when foldid is missing.
#' @param alignment Determines whether the predictions in each fold are aligned according to the fraction
#'    of progress along each fold ("fraction", default) or to use the lambda values computed on the master fit
#'    to line up the data on each fold ("lambda"). See \link[glmnet]{cv.glmnet} for details. One of
#'    c("fraction", "lambda"). Default: "fraction".
#' @param grouped Not implemented in this version. Default: FALSE.
#' @param keep If keep=TRUE, returns a list of fitted values for each fold and a corresponding list of foldids.
#'    Default: FALSE.
#' @param use_houwelingen Whether to use the houwelingen method for calculating
#' @param \dots Additional parameters to pass to the rpair_gloss or rpair_hinge function.
#'
#' @return An object of class "cv_rpair" containing the following list of values:
#'
#' @export
cv_rpair <- function(x,
                     y,
                     loss_type = c("exp", "log", "sqh", "huh"),
                     lambda=NULL,
                     nlambda=100,
                     type.measure = c("deviance", "cindex"),
                     nfolds = 5,
                     foldid = NULL,
                     alignment = c("fraction", "lambda"),
                     grouped = FALSE,
                     keep = FALSE,
                     use_houwelingen = TRUE,
                     ...
){

  loss_type=match.arg(loss_type)
  type.measure = match.arg(type.measure)
  alignment = match.arg(alignment)
  if (!is.null(lambda) && length(lambda) < 2)
    stop("Need more than one value of lambda for cv_rpair")
  if (!is.null(lambda) && alignment == "fraction") {
    warning("fraction of path alignment not available if lambda given as argument; switched to alignment=`lambda`")
    alignment = "lambda"
  }
  if(alignment == "lambda" && is.null(lambda)){
    warning("Unusual behavior has been observed when using alignment=`lambda`, particularly for logistic loss.
 Proper functionality cannot be guaranteed when selecting alignment=`lambda`, proceed at own risk.")
  }
  N = nrow(x)
  y = drop(y)
  cv.call = rpair.call = match.call(expand.dots = TRUE)
  which = match(c("type.measure", "nfolds", "foldid", "grouped",
                  "keep"), names(rpair.call), FALSE)
  if (any(which))
    rpair.call = rpair.call[-which]
  rpair.call[[1]] = as.name("rpair")

  if (is.null(foldid)){
    is_surv = ncol(y)==2
    pairs = y_to_pairs(y)
    if(identical(pairs,y)) stop("The function cv_rpair does not support pairs as an input type without user-provided
                                fold")
    if(is_surv){
      foldid_df = get_stratified_folds(y, nfolds)
    }else{
      foldid_df = get_folds(y, nfolds)
    }
    foldid = foldid_df$fold
  }else{
    nfolds = max(foldid)
    foldid_df = NULL
  }
  if (nfolds < 3)
    stop("nfolds must be bigger than 3; nfolds=10 recommended")
  cv_rpair_raw(x, y, loss_type, lambda, nlambda, type.measure,
                   nfolds, foldid, alignment, grouped, keep,
                   rpair.call, cv.call, foldid_df, use_houwelingen, ...)

}



cv_rpair_raw <- function(x, y, loss_type, lambda, nlambda, type.measure, nfolds,
                         foldid, alignment, grouped, keep, rpair.call,
                         cv.call, foldid_df, use_houwelingen, ...){

  rpair.object = switch(loss_type,
                        log = rpair_gloss(x, y, loss_type=loss_type, lambda = lambda, nlambda = nlambda, ...),
                        exp = rpair_gloss(x, y, loss_type=loss_type, lambda = lambda, nlambda = nlambda, ...),
                        sqh = rpair_hinge(x, y, loss_type=loss_type, lambda = lambda, nlambda = nlambda, ...),
                        huh = rpair_hinge(x, y, loss_type=loss_type, lambda = lambda, nlambda = nlambda, ...))
  rpair.object$call = rpair.call
  losstype = class(rpair.object)[[1]]
  type.measure = rpair_cvtype(type.measure, losstype)

  N = nrow(x)

  if(alignment=="lambda") lambda = rpair.object$lambda
  nz = sapply(predict(rpair.object, type = "nonzero"),
              length)
  seqfolds <- seq(max(foldid))
  outlist <- switch(loss_type,
                    log = lapply( seqfolds, function(i)
                      rpair_gloss(x[foldid != i,], y[foldid != i,], lambda = lambda, nlambda = nlambda, ...)),
                    exp = lapply( seqfolds, function(i)
                      rpair_gloss(x[foldid != i,], y[foldid != i,], lambda = lambda, nlambda = nlambda, ...)),
                    sqh = lapply( seqfolds, function(i)
                      rpair_hinge(x[foldid != i,], y[foldid != i,], lambda = lambda, nlambda = nlambda, ...)),
                    huh = lapply( seqfolds, function(i)
                      rpair_hinge(x[foldid != i,], y[foldid != i,], lambda = lambda, nlambda = nlambda, ...)))


  lambda = rpair.object$lambda
  class(outlist) = paste0(losstype, "list")
  predmat = rpair_buildPredmat(outlist, nlambda, lambda, x, foldid,
                         alignment)

  # calculate both deviance and concordance - type.measure will determine which one to plot
  #  deviance
  cvfl = cv_deviance(predmat, y, loss_type, foldid, delta = rpair.object$delta, use_houw = use_houwelingen)
  dev_out = cv_stats(cvfl, lambda, nz)
  lamin = with(dev_out, getopt_cv_rpair(lambda, cvm, cvsd, rpair_cvtype("deviance", losstype)))
  dev_out = c(dev_out, as.list(lamin))
  #  concordance
  cvfl = cv_concordance(predmat, y, foldid)
  conc_out = cv_stats(cvfl, lambda, nz, conc=T)
  lamin = with(conc_out, getopt_cv_rpair(lambda, cvm, cvsd, rpair_cvtype("cindex", losstype)))
  conc_out = c(conc_out, as.list(lamin))

  # assemble cross-validation results
  cvname = names(type.measure)
  names(cvname) = type.measure
  out = c(list(dev=dev_out, conc=conc_out), list(call = cv.call, name = cvname, rpair.fit = rpair.object, houw=use_houwelingen))
  if (keep)
    out = c(out, list(fit.preval = predmat, foldid = foldid, foldid_df = foldid_df))
  class(out) <- "cv_rpair"

  out

}


