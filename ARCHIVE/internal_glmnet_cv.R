# cv.glmnet
kc.cv.glmnet <- function (x, y, weights = NULL, offset = NULL, lambda = NULL, 
          type.measure = c("default", "mse", "deviance", "class", "auc", 
                           "mae", "C"), nfolds = 10, foldid = NULL, alignment = c("lambda", 
                                                                                  "fraction"), grouped = TRUE, keep = FALSE, parallel = FALSE, 
          gamma = c(0, 0.25, 0.5, 0.75, 1), relax = FALSE, trace.it = 0, 
          ...) 
{
  type.measure = match.arg(type.measure)
  alignment = match.arg(alignment)
  if (!is.null(lambda) && length(lambda) < 2) 
    stop("Need more than one value of lambda for cv.glmnet")
  if (!is.null(lambda) && alignment == "fraction") {
    warning("fraction of path alignment not available if lambda given as argument; switched to alignment=`lambda`")
    alignment = "lambda"
  }
  N = nrow(x)
  if (is.null(weights)) 
    weights = rep(1, N)
  else weights = as.double(weights)
  y = drop(y)
  cv.call = glmnet.call = match.call(expand.dots = TRUE)
  which = match(c("type.measure", "nfolds", "foldid", "grouped", 
                  "keep"), names(glmnet.call), FALSE)
  if (any(which)) 
    glmnet.call = glmnet.call[-which]
  glmnet.call[[1]] = as.name("glmnet")
  # remove trace.it code
  if (is.null(foldid)) 
    foldid = sample(rep(seq(nfolds), length = N))
  else nfolds = max(foldid)
  if (nfolds < 3) 
    stop("nfolds must be bigger than 3; nfolds=10 recommended")
  # removed relax=TRUE if-case
  kc.cv.glmnet.raw(x, y, weights, offset, lambda, type.measure, 
                nfolds, foldid, alignment, grouped, keep, parallel, trace.it, 
                glmnet.call, cv.call, ...)
}


# cv.glmnet.raw
kc.cv.glmnet.raw <- function (x, y, weights, offset, lambda, type.measure, nfolds, 
          foldid, alignment, grouped, keep, parallel, trace.it, glmnet.call, 
          cv.call, ...) 
{
  # removed trace.it code
  glmnet.object = glmnet::glmnet(x, y, weights = weights, offset = offset, 
                         lambda = lambda, trace.it = trace.it, ...)
  glmnet.object$call = glmnet.call
  subclass = class(glmnet.object)[[1]]
  type.measure = glmnet:::cvtype(type.measure, subclass)
  is.offset = glmnet.object$offset
  if (inherits(glmnet.object, "multnet") && !glmnet.object$grouped) {
    nz = predict(glmnet.object, type = "nonzero")
    nz = sapply(nz, function(x) sapply(x, length))
    nz = ceiling(apply(nz, 1, median))
  }
  else nz = sapply(predict(glmnet.object, type = "nonzero"), 
                   length)
  outlist = as.list(seq(nfolds))
  N = nrow(x)
  if (parallel) {
    outlist = foreach(i = seq(nfolds), .packages = c("glmnet")) %dopar% 
      {
        which = foldid == i
        if (length(dim(y)) > 1) 
          y_sub = y[!which, ]
        else y_sub = y[!which]
        if (is.offset) 
          offset_sub = as.matrix(offset)[!which, ]
        else offset_sub = NULL
        glmnet(x[!which, , drop = FALSE], y_sub, lambda = lambda, 
               offset = offset_sub, weights = weights[!which], 
               ...)
      }
  }
  else {
    for (i in seq(nfolds)) {
      if (trace.it) 
        cat(sprintf("Fold: %d/%d\n", i, nfolds))
      which = foldid == i
      if (is.matrix(y)) 
        y_sub = y[!which, ]
      else y_sub = y[!which]
      if (is.offset) 
        offset_sub = as.matrix(offset)[!which, ]
      else offset_sub = NULL
      outlist[[i]] = glmnet(x[!which, , drop = FALSE], 
                            y_sub, lambda = lambda, offset = offset_sub, 
                            weights = weights[!which], trace.it = trace.it, 
                            ...)
    }
  }
  lambda = glmnet.object$lambda
  class(outlist) = paste0(subclass, "list")
  predmat = buildPredmat(outlist, lambda, x, offset, foldid, 
                         alignment, y = y, weights = weights, grouped = grouped, 
                         type.measure = type.measure, family = family(glmnet.object))
  fun = paste("kc.cv", subclass, sep = ".")
  cvstuff = do.call(fun, list(predmat, y, type.measure, weights, 
                              foldid, grouped))
  grouped = cvstuff$grouped
  if ((N/nfolds < 3) && grouped) {
    warning("Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold", 
            call. = FALSE)
    grouped = FALSE
  }
  out = cvstats(cvstuff, foldid, nfolds, lambda, nz, grouped)
  cvname = names(cvstuff$type.measure)
  names(cvname) = cvstuff$type.measure
  out = c(out, list(call = cv.call, name = cvname, glmnet.fit = glmnet.object))
  if (keep) 
    out = c(out, list(fit.preval = predmat, foldid = foldid))
  lamin = with(out, getOptcv.glmnet(lambda, cvm, cvsd, cvname))
  obj = c(out, as.list(lamin))
  class(obj) = "cv.glmnet"
  obj
}

kc.cv.elnet <- function (predmat, y, type.measure, weights, foldid, grouped) 
{
  N = length(y) - apply(is.na(predmat), 2, sum)
  cvraw = switch(type.measure, mse = (y - predmat)^2, deviance = (y - 
                                                                    predmat)^2, mae = abs(y - predmat))
  list(cvraw = cvraw, weights = weights, N = N, type.measure = type.measure, 
       grouped = grouped)
}

#buildPredmat.default
# {
#   if (!is.null(offset)) {
#     is.offset = TRUE
#     offset = drop(offset)
#   }
#   else is.offset = FALSE
#   predmat = matrix(NA, nrow(x), length(lambda))
#   nfolds = max(foldid)
#   nlams = double(nfolds)
#   nlambda = length(lambda)
#   for (i in seq(nfolds)) {
#     which = foldid == i
#     fitobj = outlist[[i]]
#     if (is.offset) 
#       off_sub = offset[which]
#     preds = switch(alignment, fraction = predict(fitobj, 
#                                                  x[which, , drop = FALSE], newoffset = off_sub, ...), 
#                    lambda = predict(fitobj, x[which, , drop = FALSE], 
#                                     s = lambda, newoffset = off_sub, ...))
#     nlami = min(ncol(preds), nlambda)
#     predmat[which, seq(nlami)] = preds[, seq(nlami)]
#     if (nlami < nlambda) 
#       predmat[which, seq(from = nlami, to = nlambda)] = preds[, 
#                                                               nlami]
#   }
#   rn = rownames(x)
#   sn = paste("s", seq(0, length = nlambda), sep = "")
#   dimnames(predmat) = list(rn, sn)
#   predmat
# }