
# Change the names from '.' to '_'
# also support the case where no fold-ids are given
#   after finishing up
#   cindex - call concordance fun (already implemented)
#   loss - do in R; later consider moving to C
# add support for pairs passed to this function

cv.rpair <- function(x,
                     y,
                     lambda=NULL,
                     nlambda=100,
                     type.measure = c("default", "mse", "deviance", "auc", "mae"),
                     #type.measure = c("loss", "cindex"),
                     nfolds = 10,
                     foldid = NULL,
                     alignment = c("lambda", "fraction"),
                     grouped = FALSE,
                     keep = FALSE,
                     ...
){
  
  # how many cases?
  # two loss types: "log" and "exp"
  # three data types: ranking, survival, or pairs
  #   get working for one data type first, then expand to the others
  # two alignment types: "lambda" and "fraction"
  # four type.measure types: "mse", "deviance", "auc", "mae"
  #   possibly more - need to discuss
  
  type.measure = match.arg(type.measure)
  alignment = match.arg(alignment)
  if (!is.null(lambda) && length(lambda) < 2) 
    stop("Need more than one value of lambda for cv.rpair")
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
  cv.rpair.raw(x, y, lambda, nlambda, type.measure, 
                   nfolds, foldid, alignment, grouped, keep, 
                   rpair.call, cv.call, ...)
  
}

cv.rpair.raw <- function(x, y, lambda, nlambda, type.measure, nfolds, 
                         foldid, alignment, grouped, keep, rpair.call, 
                         cv.call, ...){

  rpair.object = rpair_gloss(x, y, lambda = lambda, nlambda = nlambda, ...)
  rpair.object$call = rpair.call
  subclass = class(rpair.object)[[1]]
  #type.measure = rpair.cvtype(type.measure, subclass)
  
  nz = sapply(predict(rpair.object, type = "nonzero"), length)
  N = nrow(x)
  
  if(alignment=="lambda") lambda = rpair.object$lambda
  
  outlist <- lapply( seq(max(foldid)), function(i) 
    rpair_gloss(x[foldid != i,], y[foldid != i,], lambda = lambda, nlambda = nlambda, ...))
  
  lambda = rpair.object$lambda
  class(outlist) = paste0(subclass, "list")
  predmat = rpair.buildPredmat(outlist, nlambda, lambda, x, foldid, 
                         alignment, y = y, grouped = grouped, 
                         type.measure = type.measure, family = class(rpair.object)[1])

  out = list(call = cv.call, rpair.fit = rpair.object)
  if (keep) 
    out = c(out, list(fit.preval = predmat, foldid = foldid))

  out
  
}

# # if fraction return df with nlambda num cols - NAs for folds that
# #   don't have values beyond some nth lambda; if lambda, will return dfs with ncols equal to length of lambda
# rpair.buildPredmat <- function(outlist, nlambda, lambda, x, foldid, alignment, ...){
#   if(alignment=="lambda") nlambda = length(lambda)
#   pred_list = list()
#   nfolds = max(foldid)
#   nlams = double(nfolds)
#   for (i in seq(nfolds)) {
#     which = foldid == i
#     fitobj = outlist[[i]]
#     predmat = matrix(NA,nrow(x[which,]),nlambda)
#     preds = switch(alignment, fraction = predict(fitobj, 
#                                                  x[which, , drop = FALSE]), 
#                    lambda = predict(fitobj, x[which, , drop = FALSE], 
#                                     s = lambda))
#     nlami = min(ncol(preds), nlambda)
#     predmat[, seq(nlami)] = preds[, seq(nlami)]
#     rn = rownames(x[which,])
#     sn = paste("s", seq(0, length = nlambda), sep = "")
#     dimnames(predmat) = list(rn, sn)
#     pred_list[[i]] = predmat
#   }
#   pred_list
# }

# I edited your implementation by deleting all .[x,], 
# becasue I need prediction for all samples to be able to calculate Houwelingen cv
rpair.buildPredmat <- function(outlist, nlambda, lambda, x, foldid, alignment, ...){
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