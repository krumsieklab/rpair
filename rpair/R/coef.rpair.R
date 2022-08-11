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
#' efit = rpair_gloss(rpair::ds1_x, rpair::ds1_y)
#' coefmat = as.matrix(coef(efit))
#'
#' @method coef rpair
#'
#' @export
coef.rpair <- function(object, s=NULL){
  predict(object, s=s, type = "coefficients")
}
