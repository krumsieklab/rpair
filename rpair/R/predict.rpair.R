#' Predict or extract values from rpair object
#'
#' Predicts fitted values given a matrix of new values. Can also extract and return coefficient and nonzero lists.
#'
#' @param object A fitted rpair object.
#' @param newx A matrix of new values. This argument is not required for type=c("coefficients", "nonzero").
#' @param s Value(s) of the penalty parameter lambda at which predictions are required. Default is the entire
#'    sequence used to create the model.
#' @param type Type of prediction desired. Type "link" gives linear predictors, type "coefficients" returns a matrix
#'    of coefficients, and type "nonzero" returns a list of the indices of the nonzero coefficients.
#'
#' @return The object returned depends on the value of \code{type}. See \code{type}.
#'
#' @examples
#' efit = rpair_gloss(rpair::ds1_x, rpair::ds1_y, standardize = F, pmax = 50, loss_type = "exp")
#' newpred = predict(efit, rpair::ds1_x, type="link")
#' nonzero = predict(efit, type = "nonzero")
#'
#' @author KC
#'
#' @method predict rpair
#'
#' @export
predict.rpair <- function(object, newx, s = NULL, type = c("link",
                                                           "coefficients", "nonzero")){

  type = match.arg(type)
  if (missing(newx)) {
    if (!match(type, c("coefficients", "nonzero"), FALSE))
      stop("You need to supply a value for 'newx'")
  }

  nbeta = object$beta

  if (!is.null(s)) {
    vnames = dimnames(nbeta)[[1]]
    dimnames(nbeta) = list(NULL, NULL)
    lambda = object$lambda
    lamlist = lambda.interp(lambda, s)
    nbeta = nbeta[, lamlist$left, drop = FALSE] %*% Matrix::Diagonal(x = lamlist$frac) +
      nbeta[, lamlist$right, drop = FALSE] %*% Matrix::Diagonal(x = 1 -
                                                          lamlist$frac)
    namess = names(s)
    if (is.null(namess))
      namess = paste0("s", seq(along = s))
    dimnames(nbeta) = list(vnames, namess)
  }

  if (type == "coefficients")
    return(nbeta)
  if (type == "nonzero")
    return(nonzeroCoef(nbeta, bystep = TRUE))
  if (inherits(newx, "sparseMatrix"))
    newx = as(newx, "dgCMatrix")
  dx = dim(newx)
  p = object$dim[1]
  if (is.null(dx))
    newx = matrix(newx, 1, byrow = TRUE)
  if (ncol(newx) != p)
    stop(paste0("The number of variables in newx must be ",
                p))

  nfit = as.matrix(newx %*% nbeta)

  nfit

}
