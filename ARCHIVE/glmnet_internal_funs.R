# glmnet:::plot.glmnet
function (x, xvar = c("norm", "lambda", "dev"), label = FALSE, 
          ...) 
{
  xvar = match.arg(xvar)
  plotCoef(x$beta, lambda = x$lambda, df = x$df, dev = x$dev.ratio, 
           label = label, xvar = xvar, ...)
}

# add number of features
# glmnet:::plotCoef
# by default use theme_minimal
function (beta, norm, lambda, df, dev, label = FALSE, xvar = c("norm", 
                                                               "lambda", "dev"), xlab = iname, ylab = "Coefficients", ...) 
{
  which = nonzeroCoef(beta)
  nwhich = length(which)
  switch(nwhich + 1, `0` = {
    warning("No plot produced since all coefficients zero")
    return()
  }, `1` = warning("1 or less nonzero coefficients; glmnet plot is not meaningful"))
  beta = as.matrix(beta[which, , drop = FALSE])
  xvar = match.arg(xvar)
  switch(xvar, norm = {
    index = if (missing(norm)) apply(abs(beta), 2, sum) else norm
    iname = "L1 Norm"
    approx.f = 1
  }, lambda = {
    index = log(lambda)
    iname = "Log Lambda"
    approx.f = 0
  }, dev = {
    index = dev
    iname = "Fraction Deviance Explained"
    approx.f = 1
  })
  dotlist = list(...)
  type = dotlist$type
  if (is.null(type)) 
    matplot(index, t(beta), lty = 1, xlab = xlab, ylab = ylab, 
            type = "l", ...)
  else matplot(index, t(beta), lty = 1, xlab = xlab, ylab = ylab, 
               ...)
  atdf = pretty(index)
  prettydf = approx(x = index, y = df, xout = atdf, rule = 2, 
                    method = "constant", f = approx.f)$y
  axis(3, at = atdf, labels = prettydf, tcl = NA)
  if (label) {
    nnz = length(which)
    xpos = max(index)
    pos = 4
    if (xvar == "lambda") {
      xpos = min(index)
      pos = 2
    }
    xpos = rep(xpos, nnz)
    ypos = beta[, ncol(beta)]
    text(xpos, ypos, paste(which), cex = 0.5, pos = pos)
  }
}

# glmnet:::coef.glmnet
function (object, s = NULL, exact = FALSE, ...) 
  predict(object, s = s, type = "coefficients", exact = exact, 
          ...)

# glmnet:::predict.glmnet
function (object, newx, s = NULL, type = c("link", "response", 
                                           "coefficients", "nonzero", "class"), exact = FALSE, newoffset, 
          ...) 
{
  type = match.arg(type)
  if (missing(newx)) {
    if (!match(type, c("coefficients", "nonzero"), FALSE)) 
      stop("You need to supply a value for 'newx'")
  }
  if (exact && (!is.null(s))) {
    lambda = object$lambda
    which = match(s, lambda, FALSE)
    if (!all(which > 0)) {
      lambda = unique(rev(sort(c(s, lambda))))
      check_dots(object, ...)
      object = update(object, lambda = lambda, ...)
    }
  }
  a0 = t(as.matrix(object$a0))
  rownames(a0) = "(Intercept)"
  nbeta = methods::rbind2(a0, object$beta)
  if (!is.null(s)) {
    vnames = dimnames(nbeta)[[1]]
    dimnames(nbeta) = list(NULL, NULL)
    lambda = object$lambda
    lamlist = lambda.interp(lambda, s)
    nbeta = nbeta[, lamlist$left, drop = FALSE] %*% Diagonal(x = lamlist$frac) + 
      nbeta[, lamlist$right, drop = FALSE] %*% Diagonal(x = 1 - 
                                                          lamlist$frac)
    namess = names(s)
    if (is.null(namess)) 
      namess = paste0("s", seq(along = s))
    dimnames(nbeta) = list(vnames, namess)
  }
  if (type == "coefficients") 
    return(nbeta)
  if (type == "nonzero") 
    return(nonzeroCoef(nbeta[-1, , drop = FALSE], bystep = TRUE))
  if (inherits(newx, "sparseMatrix")) 
    newx = as(newx, "dgCMatrix")
  dx = dim(newx)
  p = object$dim[1]
  if (is.null(dx)) 
    newx = matrix(newx, 1, byrow = TRUE)
  if (ncol(newx) != p) 
    stop(paste0("The number of variables in newx must be ", 
                p))
  nfit = as.matrix(cbind2(1, newx) %*% nbeta)
  if (object$offset) {
    if (missing(newoffset)) 
      stop("No newoffset provided for prediction, yet offset used in fit of glmnet", 
           call. = FALSE)
    if (is.matrix(newoffset) && inherits(object, "lognet") && 
        dim(newoffset)[[2]] == 2) 
      newoffset = newoffset[, 2]
    nfit = nfit + array(newoffset, dim = dim(nfit))
  }
  nfit
}

# glmnet:::pedict.glmnetfit
function (object, newx, s = NULL, type = c("link", "response", 
                                           "coefficients", "nonzero"), exact = FALSE, newoffset, ...) 
{
  type = match.arg(type)
  nfit <- NextMethod("predict")
  if (type == "response") {
    object$family$linkinv(nfit)
  }
  else {
    nfit
  }
}

# glmnet:::nonzeroCoef
function (beta, bystep = FALSE) 
{
  nr = nrow(beta)
  if (nr == 1) {
    if (bystep) 
      apply(beta, 2, function(x) if (abs(x) > 0) 
        1
        else NULL)
    else {
      if (any(abs(beta) > 0)) 
        1
      else NULL
    }
  }
  else {
    beta = abs(beta) > 0
    which = seq(nr)
    ones = rep(1, ncol(beta))
    nz = as.vector((beta %*% ones) > 0)
    which = which[nz]
    if (bystep) {
      if (length(which) > 0) {
        beta = as.matrix(beta[which, , drop = FALSE])
        nzel = function(x, which){
          if (any(x)) 
            which[x]
          else NULL
        }
        which = apply(beta, 2, nzel, which)
        if (!is.list(which)) 
          which = data.frame(which)
        which
      }
      else {
        dn = dimnames(beta)[[2]]
        which = vector("list", length(dn))
        names(which) = dn
        which
      }
    }
    else which
  }
}