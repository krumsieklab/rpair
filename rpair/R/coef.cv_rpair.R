#' Extract coefficients from a cv_rpair object
#'
#' This function is equivalent to calling predict(object, type="coefficients") on a cv_rpair object.
#'
#' @param object A fitted rpair object.
#' @param s Value(s) of the penalty parameter lambda at which predictions are required. Default is the value
#'    s="lambda.1se".
#'
#' @returns A matrix of coefficients for the provided lambda.
#'
#' @examples
#' efit = cv_rpair(ds1_x, ds1_y, nlambda=25, loss_type="exp")
#' ec = coef(efit, s="lambda.1se")
#'
#' @method coef cv_rpair
#'
#' @export
coef.cv_rpair <- function (object, s = c("lambda.1se", "lambda.min"))
{
  # determine user-selected type measure and method
  type_measure <- object$type_measure
  houw <- object$houw
  if(type_measure == "cindex" && houw){
    res_obj <- object[["h_conc"]]
  }else if(type_measure == "cindex" && !houw){
    res_obj <- object[["s_conc"]]
  }else if(type_measure == "deviance" && houw){
    res_obj <- object[["h_dev"]]
  }else if(type_measure == "deviance" && !houw){
    res_obj <- object[["s_dev"]]
  }else{
    stop("Unrecognized value(s) for type_measure and / or houw in cv_rpair object.")
  }

  # extract lambda of interest
  if (is.numeric(s))
    lambda = s
  else if (is.character(s)) {
    s = match.arg(s)
    lambda = res_obj[[s]]
  }
  else stop("Invalid form for s")

  # call coef.rpair on master fit
  coef(object$rpair.fit, s = lambda)
}
