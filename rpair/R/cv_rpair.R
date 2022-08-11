#' Cross-Validation function for rpair.
#'
#' Performs k-fold cross-validation of rpair.
#'
#' @param x Input matrix of dimension nobs x nvars; each row is an observation vector.
#' @param y Pairwise ranking analysis response variable. This function supports three types of input types:
#'    (1) continuous values, (2) survival data, and (3) ranked pairs.
#' @param loss_type Loss function to use. One of c("exp", "log" ,"sqh", "huh"). Default: "exp".
#' @param lambda Optional user-supplied lambda sequence. Default: NULL.
#' @param nlambda Number of lambda values. Default: 100.
#' @param type_measure Loss to use for cross-validation. Available loss functions are "deviance"
#'    (use loss function) and "cindex" (use concordance). Default: "deviance".
#' @param nfolds Number of folds. Default: 5.
#' @param foldid Optional vector of values between 1 and nfolds identifying which fold each observation should be in.
#'   If provided, no folds can be missing. Pairs cannot be provided as input when foldid is
#'   missing. Default: NULL.
#' @param alignment Determines whether the predictions in each fold are aligned according to the fraction
#'    of progress along each fold ("fraction", default) or to use the lambda values computed on the master fit
#'    to line up the data on each fold ("lambda"). See \link[glmnet]{cv.glmnet} for details. One of
#'    c("fraction", "lambda"). Default: "fraction".
#' @param grouped Not implemented in this version. Default: FALSE.
#' @param keep If keep=TRUE, returns a list of fitted values for each fold and a corresponding list of foldids.
#'    Default: FALSE.
#' @param use_houwelingen Whether to use the houwelingen method for calculating loss. Default: TRUE.
#' @param \dots Additional parameters to pass to the rpair function.
#'
#' @return An object of class "cv_rpair" containing the following list of values:
#' \itemize{
#'    \item{h_dev - Deviance measures calculated using houwelingen method. This is a named list
#'       containing the following:
#'       \itemize{
#'          \item{lambda - the values of lambda used in the fits.}
#'          \item{cvm - the mean cross-validated error (vector of length lambda).}
#'          \item{cvsd - estimate standard error of cvm (vector of length lambda).}
#'          \item{cvup - upper curve = cvm+cvsd.}
#'          \item{cvlo - lower curve = cvm-cvsd.}
#'          \item{nzero - number of non-zero coefficients at each lambda.}
#'          \item{lambda.min - value of lambda that produces minimum cvm.}
#'          \item{lambda.1se - largest value of lambda that error is within 1 std error of the minimum.}
#'          \item{index - a one column matrix with indices of lambda.min and lambda.1se.}
#'      }
#'    }
#'    \item{h_conc - Concordance measures calculated using houwelingen method. See h_dev.}
#'    \item{s_dev - Deviance measures calculated using standard method. See h_dev.}
#'    \item{s_conc - Concordance measures calculated using standard method. See h_dev.}
#'    \item{call - The call that produced the cv_rpair object.}
#'    \item{name - Text string indicating the loss function (for plotting purposes).}
#'    \item{type_measure - Loss to use with util functions (e.g. plot).}
#'    \item{rpair.fit - A fitted rpair object for the full data.}
#'    \item{houw - Whether to use measures calculated using houwelingen method with util functions
#'       (e.g. plot).}
#'    \item{fit.prevail - If keep=TRUE, a list of prevalidated fits for each fold. Some values
#'       can be NA if the length of lambdas calculated for that fold are shorter than the master
#'       fit.}
#'    \item{foldid - If keep=TRUE, the fold assignments used.}
#'    \item{foldid_df - Data frame consisting of columns: index, (outcome columns), and fold.}
#' }
#'
#' @export
cv_rpair <- function(x,
                     y,
                     loss_type = c("exp", "log", "sqh", "huh"),
                     lambda=NULL,
                     nlambda=100,
                     type_measure = c("deviance", "cindex"),
                     nfolds = 5,
                     foldid = NULL,
                     alignment = c("fraction", "lambda"),
                     grouped = FALSE,
                     keep = FALSE,
                     use_houwelingen = TRUE,
                     ...
){

  # verify arguments
  loss_type=match.arg(loss_type)
  type_measure = match.arg(type_measure)
  alignment = match.arg(alignment)

  # perform argument checks
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

  # cv.call - record of the function call saved to the cv_rpair object
  # rpair.call - record of the internal rpair function call made on the full dataset inside
  #   the internal cv_rpair_run function; this overwrites the returned rpair call object
  cv.call = rpair.call = match.call(expand.dots = TRUE)
  which = match(c("type_measure", "nfolds", "foldid", "grouped",
                  "keep"), names(rpair.call), FALSE)
  if (any(which))
    rpair.call = rpair.call[-which]
  rpair.call[[1]] = as.name("rpair")

  # if foldid not user-provided, generate random folds. If outcome is survival, use stratified fold
  #   function to ensure even distribution of outcomes per fold.
  y = drop(y)
  pairs = y_to_pairs(y)
  is_pairs = identical(pairs,y)
  if (is.null(foldid)){
    is_surv = is_survival(y)
    if(is_pairs) stop("The function cv_rpair does not support pairs as an input type without user-provided
                                fold ids.")
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

  cv_rpair_run(x, y, loss_type, lambda, nlambda, type_measure,
                   nfolds, foldid, alignment, grouped, keep,
                   rpair.call, cv.call, foldid_df, use_houwelingen, is_pairs, ...)

}


#' Internal function that performs cross-validation
#'
#' This function is called by cv_rpair after initial argument checks and data preparation. Contains
#' the actual code for performing cross-validation on an rpair object.
#'
#' @param x Input matrix of dimension nobs x nvars; each row is an observation vector.
#' @param y Pairwise ranking analysis response variable.
#' @param loss_type Loss function to use.
#' @param lambda Optional user-supplied lambda sequence.
#' @param nlambda Number of lambda values.
#' @param type_measure Loss to use for cross-validation. Available loss functions are "deviance"
#'    (use loss function) and "cindex" (use concordance).
#' @param nfolds Number of folds.
#' @param foldid Optional vector of values between 1 and nfolds identifying which fold each
#'    observation should be in.
#' @param alignment Determines whether the predictions in each fold are aligned according to the fraction
#'    of progress along each fold ("fraction", default) or to use the lambda values computed on the master fit
#'    to line up the data on each fold ("lambda").
#' @param gropued Not implemented in this version.
#' @param keep If keep=TRUE, returns a list of fitted values for each fold and a corresponding list
#'    of foldids.
#' @param rpair.call Record of the function call. Overwrites the call returned internally by rpair.
#' @param cv.call Record of the function call saved to cv_rpair object.
#' @param foldid_df Data frame consisting of columns: index, (outcome columns), and fold ids.
#' @param use_houwelingen Whether to use the houwelingen method for calculating loss.
#' @param is_pairs Whether the response variable is comparable pairs. This requires special
#'    treatment when assigning fold ids.
#' @param \dots Any additional arguments to pass to the rpair function.
#'
#' @return the cv_rpair object
#'
#' @noRd
cv_rpair_run <- function(x, y, loss_type, lambda, nlambda, type_measure, nfolds,
                         foldid, alignment, grouped, keep, rpair.call,
                         cv.call, foldid_df, use_houwelingen, is_pairs, ...){

  # fit model on entire dataset to obtain 'master fit' and corresponding lambda values
  rpair.object = rpair(x, y, loss_type=loss_type, lambda = lambda, nlambda = nlambda, ...)

  # calling rpair internally causes function def to be returned for call; overwrite with function
  #   call from line 66
  rpair.object$call = rpair.call
  losstype = class(rpair.object)[[1]]
  user_type_measure = type_measure
  type_measure = rpair_cvtype(type_measure, losstype)

  # perform cross-validation on each fold
  if(alignment=="lambda") lambda = rpair.object$lambda
  nz = sapply(predict(rpair.object, type = "nonzero"),
              length)
  seqfolds <- seq(max(foldid))
  if(is_pairs){
    cp_list <- lapply(seqfolds, function(i){y[y[,1] %in% which(foldid != i) & y[,2] %in% which(foldid != i),]})
    outlist <- lapply(seqfolds, function(i){x=x[foldid != i,]; if(is.null(ncol(y))){y=y[foldid != i]}else{y=cp_list[[i]]}; rpair(x=x, y=y, loss_type=loss_type, lambda = lambda, nlambda = nlambda, ...)})
  }else{
    outlist <- lapply(seqfolds, function(i){x=x[foldid != i,]; if(is.null(ncol(y))){y=y[foldid != i]}else{y=y[foldid != i,]}; rpair(x=x, y=y, loss_type=loss_type, lambda = lambda, nlambda = nlambda, ...)})
  }

  # build prediction matrix for each fold
  lambda = rpair.object$lambda
  class(outlist) = paste0(losstype, "list")
  predmat = rpair_buildPredmat(outlist, nlambda, lambda, x, foldid,
                         alignment)

  # for convenience, calculate all possible combination of method (houwelingen or standard) and metric
  #   (concordance and deviance)
  metric_list <- list(h_dev = list(houw=T, type_measure="deviance"),
                    h_conc = list(houw=T, type_measure="cindex"),
                    s_dev = list(houw=F, type_measure="deviance"),
                    s_conc = list(houw=F, type_measure="cindex"))
  out <- lapply(names(metric_list), function(x){
    type_measure <- metric_list[[x]]$type_measure
    houw <- metric_list[[x]]$houw
    if(type_measure == "deviance"){
      cvfl = cv_deviance(predmat, y, loss_type, foldid, delta = rpair.object$delta, use_houw = houw)
      out = cv_stats(cvfl, lambda, nz)
    }else{
      cvfl = cv_concordance(predmat, y, foldid, use_houw=houw)
      out = cv_stats(cvfl, lambda, nz, conc=T)
    }
    lamin = with(out, getopt_cv_rpair(lambda, cvm, cvsd, rpair_cvtype(user_type_measure, losstype)))
    out = c(out, as.list(lamin))
  })
  names(out) <- names(metric_list)

  # assemble cross-validation results
  cvname = names(type_measure)
  names(cvname) = type_measure
  out = c(out, list(call = cv.call, name = cvname, type_measure = user_type_measure, rpair.fit = rpair.object,
                    houw=use_houwelingen))
  if (keep)
    out = c(out, list(fit.preval = predmat, foldid = foldid, foldid_df = foldid_df))
  class(out) <- "cv_rpair"

  out

}


