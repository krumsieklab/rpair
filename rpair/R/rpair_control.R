#' Change fraction deviance when any upper or lower limits equal to zero
#'
#' Minimal copy of glmnet.control from the \code{glment} package. It allows internal parameters to be viewed or the
#' parameter fdev to be modified.
#'
#' @param fdev Minimum fractional change in deviance for stopping path. Default: 1e-05.
#'
#' @return A list of named values with single element: fdev.
#'
#' @noRd
rpair_control <- function (fdev = 1e-05)
{

  if (!missing(fdev))
    .Fortran("chg_fract_dev", as.double(fdev))
  value= .Fortran("get_int_parms", fdev = double(1),
                  eps = double(1), big = double(1), mnlam = integer(1),
                  devmax = double(1), pmin = double(1), exmx = double(1))
  value

}
