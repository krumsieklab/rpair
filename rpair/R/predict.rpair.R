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
#' @return The object returned depends on the value of "type".
#'
#' @examples
#' fp<- function(S){
#' time <- S[,1]
#' status <- S[,2]
#' N = length(time)
#' # for tied times
#' time[status == 0] = time[status == 0]+1e-4
#' dtimes <- time
#' dtimes[status == 0] = Inf
#' which(outer(time, dtimes, ">"), arr.ind = T)
#' }
#' # generate some random data
#' set.seed(41)
#' x = matrix(rnorm(40000),ncol = 200 )
#' S = cbind(sample(nrow(x)), rbinom(nrow(x),1,prob = 0.7))
#' # generate pairs
#' cp = fp(S)
#' efit = rpair_gloss(x, cp, standardize = F, pmax = 50, loss_type = "exp")
#' newpred = predict(efit, x, type="link")
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

# glmnet function
lambda.interp <- function (lambda, s)
{
  if (length(lambda) == 1) {
    nums = length(s)
    left = rep(1, nums)
    right = left
    sfrac = rep(1, nums)
  }
  else {
    k = length(lambda)
    sfrac <- (lambda[1] - s)/(lambda[1] - lambda[k])
    lambda <- (lambda[1] - lambda)/(lambda[1] - lambda[k])
    sfrac[sfrac < min(lambda)] <- min(lambda)
    sfrac[sfrac > max(lambda)] <- max(lambda)
    coord <- approx(lambda, seq(lambda), sfrac)$y
    left <- floor(coord)
    right <- ceiling(coord)
    sfrac = (sfrac - lambda[right])/(lambda[left] - lambda[right])
    sfrac[left == right] = 1
    sfrac[abs(lambda[left] - lambda[right]) < .Machine$double.eps] = 1
  }
  list(left = left, right = right, frac = sfrac)
}

#glmnet function
nonzeroCoef <- function (beta, bystep = FALSE)
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
        nzel = function(x, which) if (any(x))
          which[x]
        else NULL
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
