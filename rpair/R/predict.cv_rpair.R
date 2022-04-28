#' Predict or extract values from cv_rpair object
#'
#' Predicts fitted values from a cross-validated rapir object given a matrix of new values. Can also extract and
#' return coefficient and nonzero lists.
#'
#' @param object A fitted cv_rpair object.
#' @param newx A matrix of new values. This argument is not required for type=c("coefficients", "nonzero").
#' @param s Value(s) of the penalty parameter lambda at which predictions are required. Default is the value
#'    s="lambda.1se".
#' @param \dots Additional parameters to pass to the predict function.
#'
#' @return The object returned depends on the values passed via ... .
#'
#' @author KC
#'
#' @method predict cv_rpair
#'
#' @export
predict.cv_rpair <- function (object, newx, s = c("lambda.1se", "lambda.min"), ...)
{
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
  if (is.numeric(s))
    lambda = s
  else if (is.character(s)) {
    s = match.arg(s)
    lambda = res_obj[[s]]
    names(lambda) = s
  }
  else stop("Invalid form for s")
  predict(object$rpair.fit, newx, s = lambda, ...)
}
