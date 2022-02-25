#' Plot the cross-validation curve for cv.rpair
#'
#' Plots the cross-validation curve, and upper and lower standard error curves, as a function of the lambda values
#' used. Can return a ggplot object or use base R.
#'
#' @param cvobj Fitted "cv_rpair" object.
#' @param sign.lambda Either plot against log(lambda) or its negative. Default: 1.
#' @param ggplot Whether to plot and return a ggplot object. Default: F.
#'
#' @return If ggplot=T, a ggplot object.
#'
#' @author KC
#'
#' @method plot cv_rpair
#'
#' @import ggplot2
#'
#' @export
plot.cv_rpair <- function (cvobj, sign.lambda = 1, ggplot=T)
{
  if(cvobj$houw){
    plot(cvobj$cvm[-1], type="l", col="red", ylim=c(0.5,1), ylab="Houwelingen Loss")
    lines(cvobj$cvup[-1], lty=2)
    lines(cvobj$cvlo[-1], lty=2)
    return()
  }

  if(ggplot){
    plot_df <- data.frame(cvm = cvobj$cvm, cvup = cvobj$cvup, cvlo = cvobj$cvlo, nzero = cvobj$nzero)
    plot_df <- cbind(plot_df, log_lambda = sign.lambda*log(cvobj$lambda))
    plot_df <- plot_df[!is.na(plot_df$log_lambda),]
    text_y <- max(plot_df$cvup)+max(plot_df$cvup)*0.2

    g <- ggplot(plot_df, aes(x=log_lambda,y=cvm)) +
      geom_errorbar(aes(ymin=cvlo, ymax=cvup)) +
      geom_point(color="red") +
      geom_vline(xintercept = log(cvobj$lambda.min), linetype="dashed") +
      geom_vline(xintercept = log(cvobj$lambda.1se), linetype="dashed") +
      xlab("Log(lambda)") +
      ylab(cvobj$name) +
      theme_minimal() +
      annotate("text", x=plot_df$log_lambda, y=text_y, size=5,label=plot_df$nzero)

    g
  }else{
    xlab = expression(Log(lambda))
    if (sign.lambda < 0)
      xlab = paste("-", xlab, sep = "")
    plot.args = list(x = sign.lambda * log(cvobj$lambda), y = cvobj$cvm,
                     ylim = range(cvobj$cvup, cvobj$cvlo), xlab = xlab, ylab = cvobj$name,
                     type = "n")
    do.call("plot", plot.args)
    error_bars(sign.lambda * log(cvobj$lambda), cvobj$cvup, cvobj$cvlo,
               width = 0.01, col = "darkgrey")
    points(sign.lambda * log(cvobj$lambda), cvobj$cvm, pch = 20,
           col = "red")
    axis(side = 3, at = sign.lambda * log(cvobj$lambda), labels = paste(cvobj$nz),
         tick = FALSE, line = 0)
    abline(v = sign.lambda * log(cvobj$lambda.min), lty = 3)
    abline(v = sign.lambda * log(cvobj$lambda.1se), lty = 3)
    invisible()
  }
}



error_bars <- function (x, upper, lower, width = 0.02, ...)
{
  xlim <- range(x)
  barw <- diff(xlim) * width
  segments(x, upper, x, lower, ...)
  segments(x - barw, upper, x + barw, upper, ...)
  segments(x - barw, lower, x + barw, lower, ...)
  range(upper, lower)
}
