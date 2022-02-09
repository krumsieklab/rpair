# type 'response' and 'class' not currenlty included
# removed s and exact functionality
# removed ...
predict.rpair <- function(object, newx, s = NULL, type = c("link", 
                                                  "coefficients", "nonzero")){
  
  print("Hi! from predict.rpair!")
  
  type = match.arg(type)
  if (missing(newx)) {
    if (!match(type, c("coefficients", "nonzero"), FALSE)) 
      stop("You need to supply a value for 'newx'")
  }
  
  # KC: remove exact functionality for now
  # if (exact && (!is.null(s))) {
  #   lambda = object$lambda
  #   which = match(s, lambda, FALSE)
  #   if (!all(which > 0)) {
  #     lambda = unique(rev(sort(c(s, lambda))))
  #     check_dots(object, ...)
  #     object = update(object, lambda = lambda, ...)
  #   }
  # }
  
  nbeta = object$beta
  
  #KC: remove s for now
  if (!is.null(s)) {
    vnames = dimnames(nbeta)[[1]]
    dimnames(nbeta) = list(NULL, NULL)
    lambda = object$lambda
    lamlist = glmnet:::lambda.interp(lambda, s)
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
    return(glmnet:::nonzeroCoef(nbeta, bystep = TRUE))
  if (inherits(newx, "sparseMatrix")) 
    newx = as(newx, "dgCMatrix")
  dx = dim(newx)
  p = object$dim[1]
  if (is.null(dx)) 
    newx = matrix(newx, 1, byrow = TRUE)
  if (ncol(newx) != p) 
    stop(paste0("The number of variables in newx must be ", 
                p))
  #nfit = as.matrix(cbind2(1, newx) %*% nbeta)
  nfit = as.matrix(newx %*% nbeta)
  
  # KC: remove offset functionality for now
  # if (object$offset) {
  #   if (missing(newoffset)) 
  #     stop("No newoffset provided for prediction, yet offset used in fit of glmnet", 
  #          call. = FALSE)
  #   if (is.matrix(newoffset) && inherits(object, "lognet") && 
  #       dim(newoffset)[[2]] == 2) 
  #     newoffset = newoffset[, 2]
  #   nfit = nfit + array(newoffset, dim = dim(nfit))
  # }
  
  nfit
  
}


coef.rpair <- function(object){
  # removed s and exact functionality
  # removed ...
  
  print("Hi from coef.rpair!")
  predict(object, type = "coefficients")
}

#' Plot Coefficients from a "rpair" object
#' 
#' Produces a coefficient profile plot of the coefficient paths for a fitted "rpair" object. Supports "norm",
#' "lambda" and "dev" plot types. Plots are created as ggplot objects.
#' 
#' @param x
#' @param xvar
#' @param label
#' @param legend
#' 
#' @returns 
#' 
plot.rpair <- function (x, xvar = c("norm", "lambda", "dev"), label = FALSE, legend = FALSE) 
{
  # KC: took out ... for now
  
  print("Hi from plot.rpair!")
  xvar = match.arg(xvar)
  print(xvar)
  plotcoef(x$beta, lambda = x$lambda, df = x$df, dev = x$dev.ratio, 
           label = label, legend = legend, xvar = xvar)
}


plotcoef <- function (beta, norm, lambda, df, dev, label = FALSE, legend = FALSE, xvar = c("norm", 
                                                               "lambda", "dev"), xlab = iname, ylab = "Coefficients") 
{
  # KC: took out ... for now
  
  print("Hi from plotcoef!")
  
  which = glmnet:::nonzeroCoef(beta)
  nwhich = length(which)
  switch(nwhich + 1, `0` = {
    warning("No plot produced since all coefficients zero")
    return()
  }, `1` = warning("1 or less nonzero coefficients; glmnet plot is not meaningful"))
  beta = as.matrix(beta[which, , drop = FALSE])
  xvar = match.arg(xvar)
  
  print(xvar)
  
  switch(xvar, norm = {
    #index = if (missing(norm)) apply(abs(beta), 2, sum) else norm
    index = apply(abs(beta), 2, sum)
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

  names(index) = colnames(beta)
  plot_df <- as.data.frame(beta) %>% 
    tibble::rownames_to_column("variable") %>% 
    pivot_longer(-variable) %>% 
    mutate(index=index[name])
  
  g <- plot_df %>% ggplot(aes(x=index,y=value,col=variable)) + 
    geom_line() +
    xlab(iname) +
    ylab(ylab) +
    theme_minimal()
  
  if(!legend){
    g <- g + theme(legend.position = "none")
  }
  if(label){
    g <- g + geom_label_repel(data=~subset(.x,index==min(index)),
                              aes(label=variable),nudge_x=-0.5)
  }
  
  g
  
  # type = NULL
  # if (is.null(type)) 
  #   matplot(index, t(beta), lty = 1, xlab = xlab, ylab = ylab, 
  #           type = "l")
  # else matplot(index, t(beta), lty = 1, xlab = xlab, ylab = ylab)
  # atdf = pretty(index)
  # prettydf = approx(x = index, y = df, xout = atdf, rule = 2, 
  #                   method = "constant", f = approx.f)$y
  # axis(3, at = atdf, labels = prettydf, tcl = NA)
  # if (label) {
  #   nnz = length(which)
  #   xpos = max(index)
  #   pos = 4
  #   if (xvar == "lambda") {
  #     xpos = min(index)
  #     pos = 2
  #   }
  #   xpos = rep(xpos, nnz)
  #   ypos = beta[, ncol(beta)]
  #   text(xpos, ypos, paste(which), cex = 0.5, pos = pos)
  #}
}


  