#' Fit GLM or SVM for pairwise survival ranking models
#'
#' Fit a generalized linear model or SVM model for pairwise survival ranking models. GLM models use either logistic
#' loss (concordance regression) or exponential loss. SVM models use squared hinge loss or huberized hinge loss.
#' Uses a modified implementation of \link[glmnet]{glmnet} for logistic and exponential loss and a modified
#' implementation of \link[gcdnet]{gcdnet} for squared and huberized hinge loss. Refer to glmnet and gcdnet
#' documentation for further details.
#'
#' @param x Input matrix of dimension nobs x nvars; each row is an observation vector.
#' @param y Pairwise ranking analysis response variable. This function supports three types of input types:
#'    (1) continuous values, (2) survival data, and (3) ranked pairs.
#' @param loss_type Loss function to use. Equivalent to "family" in glmnet package. One of c("log", "exp").
#' @param alpha For GLM models only. The elasticnet mixing parameter. See \link[glmnet]{glmnet} for details.
#'    Default: 1.0 (lasso penalty).
#' @param nlambda Number of lambda values to evaluate. Default: 100.
#' @param lambda.min.ratio For SVM models, this is equivalent to lambda.factor. Smallest value for lambda as a
#'   fraction of lambda.max, the (data derived) entry value. See \link[glmnet]{glmnet} or \link[gcdnet]{gcdnet} for
#'   details. Default: ifelse(nobs0 <= nvars, 1e-2, 1e-4)
#' @param lambda A user supplied lambda sequence. Overrides the typical usage in which a lambda sequence is computed
#'    using nlambda and lambda.min.ratio. Provide a decreasing sequence of lambda values with at least 2 entries.
#'    Default: NULL.
#' @param lambda2 For SVM models only. Regularization parameter lambda2 for the quadratic penalty of the coefficients.
#'    Default: 0.
#' @param standardize Logical flag for x variable standardization. Default: FALSE.
#' @param thresh For SVM models, this is equivalent to eps. Convergence threshold for coordinate descent.
#'    GLM Default: 1e-7, SVM Default: 1e-6.
#' @param dfmax Limit the maximum number of variables in the model. Default: nvars+1.
#' @param pmax Limit the maximum number of variables that can be nonzero. Default: min(dfmax*2+20, nvars).
#' @param penalty.factor For GLM models only. Vector of penalty factors to apply to each coefficient.
#'    Default: rep(1, nvars).
#' @param pf For SVM models only. L1 penalty factor of length p used for adaptive LASSO or adaptive elastic net. See
#'    \link[gcdnet]{gcdnet} for details. Default: rep(1, nvars).
#' @param pf2 For SVM models only. L2 penalty factor of length p used for adaptive LASSO or adaptive elastic net. See
#'    \link[gcdnet]{gcdnet} for details. Default: rep(1, nvars).
#' @param maxit Maximum number of passes over the data for all lambda values. Default: 100000.
#' @param type.logistic For GLM models only. One of c("Newton", "modified.Newton"). If "Newton" then the exact
#'    hessian is used, while "modified.Newton" uses an upper-bound on the hessian and can be faster.
#'    Default: "Newton".
#' @param delta For SVM models only. The parameter delta in the HHSVM model. Must be greater than 0. Default: 3.
#'
#' @return An object with S3 class "rpair", "*", where "*" is "lognet", "fishnet", "phuhnet". or "psqhnet". Contains
#' the following attributes:
#'
#' @examples
#'fp<- function(S){
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
#' efit = rpair(x, cp, standardize = F, pmax = 50, loss_type = "exp")
#' lfit = rpair(x, cp, standardize = F, pmax = 50, loss_type = "sqh")
#'
#' @author mubu, KC
#'
#' @export
rpair <- function(x,
                  y,
                  loss_type=c("exp", "log", "sqh", "huh"),
                  alpha,
                  nlambda,
                  lambda.min.ratio,
                  lambda,
                  lambda2,
                  standardize,
                  thresh,
                  dfmax,
                  pmax,
                  penalty.factor,
                  pf1,
                  pf2,
                  maxit,
                  type.logistic,
                  delta
){

  loss_type = match.arg(loss_type)
  args <- as.list(match.call())[-1]

  fit <- switch(loss_type,
                exp = do.call(rpair_gloss, args),
                log = do.call(rpair_gloss, args),
                sqh = do.call(rpair_hinge, args),
                huh = do.call(rpair_hinge, args))
  fit

}