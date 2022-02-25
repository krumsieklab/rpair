#' Extract coefficients from an rpair object
#'
#' This function is equivalent to calling predict(object, type="coefficients").
#'
#' @param object A fitted rpair object.
#' @param s Value(s) of the penalty parameter lambda at which predictions are required. Default is the entire
#'    sequence used to create the model.
#'
#' @returns A matrix of coefficients.
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
#' coefmat = coef(efit)
#'
#' @method coef rpair
#'
#' @export
coef.rpair <- function(object, s=NULL){
  predict(object, s=s, type = "coefficients")
}
