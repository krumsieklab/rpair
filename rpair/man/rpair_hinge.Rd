% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rpair_hinge.R
\name{rpair_hinge}
\alias{rpair_hinge}
\title{Fits regularization paths for large margin classifiers}
\usage{
rpair_hinge(
  x,
  y,
  loss_type = c("sqh", "huh"),
  nlambda = 100,
  lambda.factor = ifelse(nobs <= nvars, 0.01, 1e-04),
  lambda = NULL,
  lambda2 = 0,
  pf = rep(1, nvars),
  pf2 = rep(1, nvars),
  dfmax = nvars + 1,
  pmax = min(ncol(x), nrow(x)),
  standardize = FALSE,
  eps = 1e-06,
  maxit = 1e+05,
  delta = 3
)
}
\arguments{
\item{x}{Input matrix of dimension nobs x nvars; each row is an observation vector.}

\item{y}{Pairwise ranking analysis response variable. This function supports three types of input types:
(1) continuous values, (2) survival data, and (3) ranked pairs.}

\item{loss_type}{Loss function to use. Equivalent to "method" in gcdnet package. One of c("sqh", "huh").}

\item{nlambda}{Number of lambda values to evaluate. Default: 100.}

\item{lambda.factor}{The factor for getting the minimal lambda in lambda sequence, where
min(lambda) = lambda.factor * max(lambda). max(lambda) is the smallest value of lambda for which all
coefficients are zero. See \link[gcdnet]{gcdnet} for details. Default: ifelse(nobs <= nvars, 0.01, 1e-04).}

\item{lambda}{A user supplied lambda sequence. Overrides the typical usage in which a lambda sequence is computed
using nlambda and lambda.min.ratio. Provide a decreasing sequence of lambda values with at least 2 entries.
Default: NULL.}

\item{lambda2}{Regularization parameter lambda2 for the quadratic penalty of the coefficients. Default: 0.}

\item{pf}{L1 penalty factor of length p used for adaptive LASSO or adaptive elastic net. See
\link[gcdnet]{gcdnet} for details. Default: rep(1, nvars).}

\item{pf2}{L2 penalty factor of length p used for adaptive LASSO or adaptive elastic net. See
\link[gcdnet]{gcdnet} for details. Default: rep(1, nvars).}

\item{dfmax}{Limit the maximum number of variables in the model. Default: nvars+1.}

\item{pmax}{Limit the maximum number of variables that can be nonzero. Default: min(ncol(x), nrow(x)).}

\item{standardize}{Logical flag for x variable standardization. Default: FALSE.}

\item{eps}{Convergence threshold for coordinate descent. Default: 1e-06.}

\item{maxit}{Maximum number of passes over the data for all lambda values. Default: 1e+05}

\item{delta}{The parameter delta in the HHSVM model. Must be greater than 0. Default: 3.}
}
\value{
An object with S3 class \code{"rpair"}, "*", where "*" is \code{"psqhnet"} or \code{"phuhnet"}. Contains
the following attributes:
   \item{beta}{a nvars x length(lambda) matrix of coefficients, stored in sparse column format}
   \item{df}{the number of nonzero coefficients for each value of lambda}
   \item{dim}{dimension of coefficient matrix}
   \item{lambda}{The actual sequence of lambda values used}
   \item{npasses}{total passes over the data summed over all lambda values}
   \item{jerr}{error flag, for warnings and errors (largely for internal debugging)}
   \item{call}{the call that produced the object}
   \item{loss}{the loss function used}
   \item{nobs}{the number of observations}
}
\description{
Fit a regularization path for large margin classifiers at a sequence of regularization parameters lambda. Uses
a modified implementation of \link[gcdnet]{gcdnet}. Only supports squared hinge loss and huberized square hinge
loss. Refer to gcdnet documentation for further details.
}
\examples{
sfit = rpair_hinge(ds1_x, ds1_y, loss_type = "sqh")
hfit = rpair_hinge(ds1_x, ds1_y, loss_type = "huh")

}
\author{
mubu, KC
}
