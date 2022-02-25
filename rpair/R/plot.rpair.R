#' Plot coefficients from a "rpair" object
#'
#' Produces a coefficient profile plot of the coefficient paths for a fitted "rpair" object. Supports "norm",
#' "lambda" and "dev" plot types. Plots are created as ggplot objects.
#'
#' @param x An "rpair" object.
#' @param xvar Plot type. One of c("norm", "lambda", or "dev").
#' @param label Whether to include variable labels. Default: FALSE.
#' @param legend Whether to include a legend. Default: FALSE.
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
#' plot(efit, xvar="norm")
#' plot(efit, xvar="lambda", label=T)
#'
#' @author KC
#'
#' @returns ggplot object
#'
#' @method plot rpair
#'
#' @export
plot.rpair <- function (x, xvar = c("norm", "lambda", "dev"), legend = FALSE)
{

  xvar = match.arg(xvar)
  plotcoef(x$beta, lambda = x$lambda, df = x$df, dev = x$dev.ratio,
           legend = legend, xvar = xvar)
}


#' Plot coefficient internal function
#'
#' This function is called by plot.rpair. It creates the rpair coefficient plot (currently separate implementation
#' in case we add additional plotting options).
#'
#' @param beta Matrix of coefficients.
#' @param lambda The acutal sequence of lambda values used for an rpair object.
#' @param df The number of nonzero coefficinets for each value of lambda.
#' @param dev The fraction of (null) deviance explained. dev.ratio from an rpair object.
#' @param label Whether to include feature labels.
#' @param legend Whether to include a legend.
#' @param xvar Plot type. One of c("norm", "lambda", "dev").
#' @param xlab Label for x-axis.
#' @param ylab Label for y-axis.
#'
#' @returns ggplot object
#'
#' @import ggplot2
#'
#' @noRd
plotcoef <- function (beta, lambda, df, dev, legend, xvar = c("norm", "lambda", "dev"),
                      xlab = iname, ylab = "Coefficients")
{

  which = nonzeroCoef(beta)
  nwhich = length(which)
  switch(nwhich + 1, `0` = {
    warning("No plot produced since all coefficients zero")
    return()
  }, `1` = warning("1 or less nonzero coefficients; rpair plot is not meaningful"))
  beta = as.matrix(beta[which, , drop = FALSE])
  xvar = match.arg(xvar)

  switch(xvar, norm = {
    index = apply(abs(beta), 2, sum)
    iname = "L1 Norm"
    approx.f = 1
  }, lambda = {
    index = log(lambda)
    iname = "Log Lambda"
    approx.f = 0
  }, dev = {
    index = dev
    iname = "Fraction Deviance Explained"
    approx.f = 1
  })

  names(index) = colnames(beta)
  plot_df <- as.data.frame(beta)
  plot_df <- cbind(variable=rownames(plot_df), plot_df)
  plot_df <- reshape(plot_df, direction = "long",
                     v.names = "value",
                     varying = 2:ncol(plot_df),
                     times = colnames(plot_df)[2:ncol(plot_df)],
                     timevar = "name",
                     idvar = "variable")
  plot_df <- cbind(index=index[plot_df$name], plot_df)
  text_y = max(plot_df$value)+max(plot_df$value)*0.3

  g <- ggplot(plot_df, aes(x=index,y=value,col=variable)) +
    geom_line() +
    xlab(iname) +
    ylab(ylab) +
    theme_minimal() +
    annotate("text", x=index, y=text_y, label=df)

  if(!legend){
    g <- g + theme(legend.position = "none")
  }


  g

}

# glmnet function
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
