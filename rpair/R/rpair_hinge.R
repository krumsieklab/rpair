#' Fits regularization paths for large margin classifiers
#'
#' Fit a regularization path for large margin classifiers at a sequence of regularization parameters lambda. Uses
#' a modified implementation of \link[gcdnet]{gcdnet}. Only supports squared hinge loss and huberized square hinge
#' loss. Refer to gcdnet documentation for further details.
#'
#' @param x Input matrix of dimension nobs x nvars; each row is an observation vector.
#' @param y Pairwise ranking analysis response variable. This function supports three types of input types:
#'    (1) continuous values, (2) survival data, and (3) ranked pairs.
#' @param loss_type Loss function to use. Equivalent to "method" in gcdnet package. One of c("sqh", "huh").
#' @param nlambda Number of lambda values to evaluate. Default: 100.
#' @param lambda.factor The factor for getting the minimal lambda in lambda sequence, where
#'    min(lambda) = lambda.factor * max(lambda). max(lambda) is the smallest value of lambda for which all
#'    coefficients are zero. See \link[gcdnet]{gcdnet} for details. Default: ifelse(nobs <= nvars, 0.01, 1e-04).
#' @param lambda A user supplied lambda sequence. Overrides the typical usage in which a lambda sequence is computed
#'    using nlambda and lambda.min.ratio. Provide a decreasing sequence of lambda values with at least 2 entries.
#'    Default: NULL.
#' @param lambda2 Regularization parameter lambda2 for the quadratic penalty of the coefficients. Default: 0.
#' @param pf L1 penalty factor of length p used for adaptive LASSO or adaptive elastic net. See
#'    \link[gcdnet]{gcdnet} for details. Default: rep(1, nvars).
#' @param pf2 L2 penalty factor of length p used for adaptive LASSO or adaptive elastic net. See
#'    \link[gcdnet]{gcdnet} for details. Default: rep(1, nvars).
#' @param dfmax Limit the maximum number of variables in the model. Default: nvars+1.
#' @param pmax Limit the maximum number of variables that can be nonzero. Default: min(ncol(x), nrow(x)).
#' @param standardize Logical flag for x variable standardization. Default: FALSE.
#' @param eps Convergence threshold for coordinate descent. Default: 1e-06.
#' @param maxit Maximum number of passes over the data for all lambda values. Default: 1e+05
#' @param delta The parameter delta in the HHSVM model. Must be greater than 0. Default: 3.
#'
#' @return An object with S3 class \code{"rpair"}, "*", where "*" is \code{"psqhnet"} or \code{"phuhnet"}. Contains
#' the following attributes:
#'    \item{beta}{a nvars x length(lambda) matrix of coefficients, stored in sparse column format}
#'    \item{df}{the number of nonzero coefficients for each value of lambda}
#'    \item{dim}{dimension of coefficient matrix}
#'    \item{lambda}{The actual sequence of lambda values used}
#'    \item{npasses}{total passes over the data summed over all lambda values}
#'    \item{jerr}{error flag, for warnings and errors (largely for internal debugging)}
#'    \item{call}{the call that produced the object}
#'    \item{loss}{the loss function used}
#'    \item{nobs}{the number of observations}
#'
#' @examples
#' sfit = rpair_hinge(surv_x, surv_cp, standardize = F, pmax = 50, loss_type = "sqh")
#' hfit = rpair_hinge(surv_x, surv_cp, standardize = F, pmax = 50, loss_type = "huh")
#'
#' @author mubu, KC
#'
#' @export
rpair_hinge <- function(x,
                        y,
                        loss_type = c("sqh", "huh"),
                        nlambda = 100,
                        lambda.factor = ifelse(nobs <= nvars, 0.01, 1e-04),
                        lambda = NULL,
                        lambda2 = 0,
                        pf = rep(1, nvars),
                        pf2 = rep(1,nvars),
                        dfmax = nvars + 1,
                        pmax= min(ncol(x), nrow(x)),
                        standardize = FALSE,
                        eps = 1e-06,
                        maxit = 1e+05,
                        delta = 3
  ){

    # return call with fitted object
    this.call <- match.call()

    ### Prepare input and outcome variables
    x <- as.matrix(x)

    nobs <- as.integer(nrow(x))
    nvars <- as.integer(ncol(x))
    vnames <- colnames(x)
    if (is.null(vnames)) vnames <- paste("V", seq(nvars), sep = "")

    # generate comparable pairs
    is_surv = ncol(y)==2
    cp = y_to_pairs(y)

    ### Prepare all the generic arguments, then hand off to loss_type functions
    ## unmodified parameters
    loss_type <- match.arg(loss_type)
    maxit <- as.integer(maxit)
    lam2 <- as.double(lambda2)
    pf <- as.double(pf)
    pf2 <- as.double(pf2)
    isd <- as.integer(standardize)
    eps <- as.double(eps)
    dfmax <- as.integer(dfmax)

    # for now, not supporting user-provided exclude
    jd <- as.integer(0)

    ## parameters with checks / conditions
    if (length(pf) != nvars)
      stop("The size of L1 penalty factor must be same as the number of input variables")
    if (length(pf2) != nvars)
      stop("The size of L2 penalty factor must be same as the number of input variables")
    if (lambda2 < 0)
      stop("lambda2 must be non-negative")

    nlam <- as.integer(nlambda)
    if (is.null(lambda)) {
      if (lambda.factor >= 1)
        stop("lambda.factor should be less than 1")
      flmin <- as.double(lambda.factor)
      ulam <- double(1)
    }
    else {
      flmin <- as.double(1)
      if (any(lambda < 0))
        stop("lambdas should be non-negative")
      ulam <- as.double(rev(sort(lambda)))
      nlam <- as.integer(length(lambda))
    }

    # use either number of features or number of non-censored data points
    if(is_surv) pmax = min(nvars, sum(y[,2]))
    pmax=as.integer(pmax)

    ### Fit model
    fit <- switch(loss_type,
                  huh = phuhnetfit(x, cp[,2:1], nlam, flmin, ulam, isd, eps, dfmax, pmax, jd, pf, pf2, maxit, lam2, delta, nobs, nvars, vnames),
                  sqh = psqhnetfit(x, cp[,2:1], nlam, flmin, ulam, isd, eps, dfmax, pmax, jd, pf, pf2, maxit, lam2, nobs, nvars, vnames)
    )

    if (is.null(lambda))  fit$lambda <- fix_lam(fit$lambda)
    fit$call <- this.call
    fit$loss=loss_type
    fit$nobs=nobs

    class(fit) <- c(class(fit), "rpair")
    fit

  }

#' Huberized squared hinge loss Function
#'
#' Helper function to call the Fortran implemented huberized squared hinge loss algorithm.
#'
#' @param x Input matrix of dimension nobs x nvars; each row is an observation vector.
#' @param cp Survival data outcome as comparable pairs.
#' @param nlam Number of lambda values to evaluate.
#' @param flmin If lambda sequence not provided, this is the lambda.min.ratio. If it is, then it is equal to 1.
#' @param ulam If provided, this is the user lambda sequence. If not, this is 1.
#' @param isd Integer value for rpair_gloss 'standardize' argument.
#' @param eps Convergence threshold for coordinate descent.
#' @param dfmax Limit the maximum number of variables in the model.
#' @param pmax Limit the maximum number of variables that can be nonzero.
#' @param jd A vector of features to be excluded.
#' @param pf L1 penalty factor of length p used for adaptive LASSO or adaptive elastic net.
#' @param pf2 L2 penalty factor of length p used for adaptive LASSO or adaptive elastic net.
#' @param maxit Maximum number of passes over the data for all lambda values.
#' @param lam2 Regularization parameter lambda2 for the quadratic penalty of the coefficients.
#' @param delta The parameter delta in the HHSVM model.
#' @param nobs Number of samples (first dimension of x matrix).
#' @param nvars Number of features (second dimension of x matrix).
#' @param vnames Column names of the input matrix.
#' @param rk Half the number of pairs.
#'
#' @author mubu, KC
#'
#' @noRd
phuhnetfit <- function( x, cp, nlam, flmin, ulam, isd, eps, dfmax, pmax, jd, pf, pf2, maxit, lam2,
                     delta, nobs, nvars, vnames, rk = nrow(cp)%/%2
){

  # modify cp accordingly
  cp = rbind( cp[1:rk,2:1], cp[-(1:rk),])
  y = c(rep(1,rk),rep(-1, nrow(cp)-rk))
  y = cbind(c0= -y, c1 =y)
  # ---

  if (delta < 0)
    stop("delta must be non-negative")
  delta <- as.double(delta)
  fit <- .Fortran("phuhnet", delta, lam2, nobs, nvars,
                  as.double(x), as.double(y), jd, pf, pf2, dfmax, pmax,
                  nlam, flmin, ulam, eps, isd, maxit, nalam = integer(1),
                  b0 = double(nlam), beta = double(pmax * nlam), ibeta = integer(pmax),
                  nbeta = integer(nlam), alam = double(nlam), npass = integer(1),
                  as.integer(nrow(cp)), as.integer(cp),
                  jerr = integer(1))
  class(fit) = c("gcdnetfit",class(fit))

  errmsg <- jerr(fit, maxit, pmax, "huh")
  switch(paste(errmsg$n), `1` = stop(errmsg$msg, call. = FALSE), `-1` = print(errmsg$msg, call. = FALSE))

  outlist <- getcoef(fit,nvars ,pmax, vnames)
  outlist <- c(outlist, list(delta=delta, npasses = fit$npass, jerr = fit$jerr))
  class(outlist) <- c("phuhnet")
  outlist
}

#' Squared hinge loss Function
#'
#' Helper function to call the Fortran implemented squared hinge loss algorithm.
#'
#' @param x Input matrix of dimension nobs x nvars; each row is an observation vector.
#' @param cp Survival data outcome as comparable pairs.
#' @param nlam Number of lambda values to evaluate.
#' @param flmin If lambda sequence not provided, this is the lambda.min.ratio. If it is, then it is equal to 1.
#' @param ulam If provided, this is the user lambda sequence. If not, this is 1.
#' @param isd Integer value for rpair_gloss 'standardize' argument.
#' @param eps Convergence threshold for coordinate descent.
#' @param dfmax Limit the maximum number of variables in the model.
#' @param pmax Limit the maximum number of variables that can be nonzero.
#' @param jd A vector of features to be excluded.
#' @param pf L1 penalty factor of length p used for adaptive LASSO or adaptive elastic net.
#' @param pf2 L2 penalty factor of length p used for adaptive LASSO or adaptive elastic net.
#' @param maxit Maximum number of passes over the data for all lambda values.
#' @param lam2 Regularization parameter lambda2 for the quadratic penalty of the coefficients.
#' @param nobs Number of samples (first dimension of x matrix).
#' @param nvars Number of features (second dimension of x matrix).
#' @param vnames Column names of the input matrix.
#' @param rk Half the number of pairs.
#'
#' @author mubu, KC
#'
#' @noRd
psqhnetfit <- function( x, cp, nlam, flmin, ulam, isd, eps, dfmax, pmax, jd, pf, pf2, maxit, lam2,
                     nobs, nvars, vnames, rk = nrow(cp)%/%2
){

  # modify cp accordingly
  cp = rbind( cp[1:rk,2:1], cp[-(1:rk),])
  y = c(rep(1,rk),rep(-1, nrow(cp)-rk))
  y = cbind(c0= -y, c1 =y)
  # ---

  fit <- .Fortran("psqhnet", lam2, nobs, nvars, as.double(x),
                  as.double(y), jd, pf, pf2, dfmax, pmax, nlam, flmin,
                  ulam, eps, isd, maxit, nalam = integer(1), b0 = double(nlam),
                  beta = double(pmax * nlam), ibeta = integer(pmax), nbeta = integer(nlam),
                  alam = double(nlam), npass = integer(1),
                  as.integer(nrow(cp)), as.integer(cp),
                  jerr = integer(1))

  class(fit) = c("gcdnetfit",class(fit))
  errmsg <- jerr(fit, maxit, pmax, "sqh")
  switch(paste(errmsg$n), `1` = stop(errmsg$msg, call. = FALSE), `-1` = print(errmsg$msg, call. = FALSE))

  outlist <- getcoef(fit,nvars ,pmax, vnames)
  outlist <- c(outlist, list(npasses = fit$npass, jerr = fit$jerr))
  class(outlist) <- c("psqhnet")
  outlist
}

