# KC: TO DO:
# - remove survival package dependency

#' Convert response variable to comparable pairs
#'
#' Convert supported response variable input types to comparable pairs. Comparable pairs of the format: two column
#' matrix where integers represent samples. Values in the left column have a greater survival time than values in
#' the right column.
#'
#' @param y Pairwise ranking analysis response variable. This function supports three types of input types:
#'    (1) continuous values, (2) survival data, and (3) ranked pairs.
#'
#' @author mubu
#'
#' @noRd
y_to_pairs <- function(y){

  require(survival)

  # type of the outcome
  if(is.null(nrow(y))){ y = as.numeric(y); ytype = "ranking"}
  else if( inherits(y,"Surv")) ytype = "survival"
  else if( identical(colnames(y),c("time","status")) ){y = Surv(time = y[,1],event = y[,2]); ytype = "survial"}
  else if( identical(colnames(y),c("start", "stop", "status")) ){y = Surv(y[,1],y[,2], y[,3]); ytype = "survial" }
  else if( ncol(y) == 2 ){ ytype = "pairs"; return(y)}
  else stop("y is not of a recognazied type for pairwise ranking")

  UseMethod("y_to_pairs", y)

}

#' Survival y-to-pairs
#'
#' Convert survival data response variable to comparable pairs.
#'
#' @param S Survival data response variable.
#'
#' @author mubu
#'
#' @noRd
y_to_pairs.Surv <- function(S){

  if(inherits(S,"Surv")){
    type <- attr(S, "type")
    if (!type %in% c("right", "counting"))
      stop("doesn't support \"", type, "\" survival data")
  }

  k = ncol(S)
  time <- S[,k-1]
  status <- S[,k]
  # for tied censored times
  time[status == 0] = time[status == 0]+1e-4
  dtimes <- time
  dtimes[status == 0] = Inf
  if (k == 2) {
    return(which(outer(time, dtimes, ">"), arr.ind = T))
  }else{
    start <- S[,1]
    return(which(outer(time, dtimes, ">") & outer(start, dtimes, "<"), arr.ind = T))
  }
}

#' Factor y-to-pairs
#'
#' Convert factor response variable to comparable pairs.
#'
#' @param y Factor response variable.
#'
#' @author mubu
#'
#' @noRd
y_to_pairs.factor <- function(y){
  if (length(unique(y)) <2)
    stop("doesn't support factor with only one level")
  q = as.numeric(y)
  return(which(outer(q, q, ">"), arr.ind = T))
}

#' Integer y-to-pairs
#'
#' Convert integer response variable to comparable pairs.
#'
#' @param y Integer response variable.
#'
#' @author mubu
#'
#' @noRd
y_to_pairs.integer <- function(q){
  if (length(unique(q)) <2)
    stop("doesn't support integer with only one value")
  return(which(outer(q, q, ">"), arr.ind = T))
}

#' Numeric y-to-pairs
#'
#' Convert continuous response variable to comparable pairs.
#'
#' @param y Continuous response variable.
#'
#' @author mubu
#'
#' @noRd
y_to_pairs.numeric <- function(q){
  if (length(unique(q)) <2)
    stop("doesn't support numeric with only one value")
  return(which(outer(q, q, ">"), arr.ind = T))
}

