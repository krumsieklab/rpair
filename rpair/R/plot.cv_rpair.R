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
  type_measure = names(cvobj$name)
  if(type_measure == "cindex"){
    res_obj <- cvobj[["conc"]]
    title <- ifelse(cvobj$houw, "Houwelingen Method", "Standard Method")
    if(ggplot){
      plot_df <- data.frame(cvm = res_obj$cvm[-1], cvup = res_obj$cvup[-1], cvlo = res_obj$cvlo[-1])
      plot_df <- cbind(index = 1:nrow(plot_df), plot_df)
      plot_df <- reshape2::melt(plot_df, id="index")

      g <- ggplot(plot_df, aes(x=index,y=value, group=variable)) +
        geom_line(aes(linetype=variable, color=variable)) +
        scale_linetype_manual(values=c("solid", "dashed", "dashed")) +
        scale_color_manual(values=c("red", "black", "black")) +
        theme_minimal() +
        ylab(cvobj$name) +
        xlab("lambda index") +
        theme(legend.position = "none") +
        ggtitle(title)

      return(g)
    }else{
      plot(res_obj$cvm[-1], type="l", col="red", ylim=c(0.5,1), ylab=cvobj$name, xlab="lambda index")
      lines(res_obj$cvup[-1], lty=2)
      lines(res_obj$cvlo[-1], lty=2)
      title(title)
      return()
    }
  }

  res_obj <- cvobj[["dev"]]
  title <- ifelse(cvobj$houw, "Houwelingen Method", "Standard Method")
  if(ggplot){
    plot_df <- data.frame(cvm = res_obj$cvm, cvup = res_obj$cvup, cvlo = res_obj$cvlo, nzero = res_obj$nzero)
    plot_df <- cbind(plot_df, log_lambda = sign.lambda*log(res_obj$lambda))
    plot_df <- plot_df[!is.na(plot_df$log_lambda),]
    text_y <- max(plot_df$cvup)+max(plot_df$cvup)*0.2

    g <- ggplot(plot_df, aes(x=log_lambda,y=cvm)) +
      geom_errorbar(aes(ymin=cvlo, ymax=cvup)) +
      geom_point(color="red") +
      geom_vline(xintercept = log(res_obj$lambda.min), linetype="dashed") +
      geom_vline(xintercept = log(res_obj$lambda.1se), linetype="dashed") +
      xlab("Log(lambda)") +
      ylab(cvobj$name) +
      theme_minimal() +
      annotate("text", x=plot_df$log_lambda, y=text_y, size=5,label=plot_df$nzero) +
      ggtitle(title)

    g
  }else{
    xlab = expression(Log(lambda))
    if (sign.lambda < 0)
      xlab = paste("-", xlab, sep = "")
    plot.args = list(x = sign.lambda * log(res_obj$lambda), y = res_obj$cvm,
                     ylim = range(res_obj$cvup, res_obj$cvlo), xlab = xlab, ylab = cvobj$name,
                     type = "n")
    do.call("plot", plot.args)
    error_bars(sign.lambda * log(res_obj$lambda), res_obj$cvup, res_obj$cvlo,
               width = 0.01, col = "darkgrey")
    points(sign.lambda * log(res_obj$lambda), res_obj$cvm, pch = 20,
           col = "red")
    axis(side = 3, at = sign.lambda * log(res_obj$lambda), labels = paste(res_obj$nzero),
         tick = FALSE, line = -1)
    abline(v = sign.lambda * log(res_obj$lambda.min), lty = 3)
    abline(v = sign.lambda * log(res_obj$lambda.1se), lty = 3)
    title(title)
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
