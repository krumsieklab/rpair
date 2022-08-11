#' Plot the cross-validation curve for cv.rpair
#'
#' Plots the cross-validation curve, and upper and lower standard error bars (deviance) /
#' curves (concordance), as a function of the lambda values. Whether deviance or
#' concordance is plotted is determined by type.measure in cv_rpair fit object. Can return
#' a ggplot object or use base R.
#'
#' @param cvobj Fitted "cv_rpair" object.
#' @param sign.lambda Either plot against log(lambda) or its negative. Default: 1.
#' @param ggplot Whether to plot and return a ggplot object. Default: T.
#'
#' @return If ggplot=T, a ggplot object. If ggplot=F, a base R plot.
#'
#' @author KC
#'
#' @method plot cv_rpair
#'
#' @import ggplot2
#'
#' @export
plot.cv_rpair <- function (cvobj, sign.lambda = 1, ggplot=T, unique_nz = F, nz_size = 3)
{
  ## plot concordance ----
  if(cvobj$type_measure == "cindex"){
    if(cvobj$houw){
      res_obj <- cvobj[["h_conc"]]
    }else{
      res_obj <- cvobj[["s_conc"]]
    }
    if(ggplot){ # ggplot
      # extract mean and standard error values
      plot_df <- data.frame(cvm = res_obj$cvm[-1], cvup = res_obj$cvup[-1], cvlo = res_obj$cvlo[-1])
      plot_df <- cbind(index = 1:nrow(plot_df), plot_df)
      plot_df <- reshape2::melt(plot_df, id="index")

      # plot curves
      g <- ggplot(plot_df, aes(x=index,y=value, group=variable)) +
        geom_line(aes(linetype=variable, color=variable)) +
        scale_linetype_manual(values=c("solid", "dashed", "dashed")) +
        scale_color_manual(values=c("red", "black", "black")) +
        theme_minimal() +
        ylab(cvobj$name) +
        xlab("lambda index") +
        theme(legend.position = "none")

      return(g)
    }else{ # base R plot
      plot(res_obj$cvm[-1], type="l", col="red", ylim=range(res_obj$cvlo, res_obj$cvup), ylab=cvobj$name, xlab="lambda index")
      lines(res_obj$cvup[-1], lty=2)
      lines(res_obj$cvlo[-1], lty=2)
      return()
    }
  }

  ## plot deviance ----
  if(cvobj$houw){
    res_obj <- cvobj[["h_dev"]]
  }else{
    res_obj <- cvobj[["s_dev"]]
  }

  if(ggplot){ # ggplot
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
      annotate("text", x=plot_df$log_lambda, y=text_y, size=nz_size,label=plot_df$nzero)

    g

  }else{ # base R plot (copied from glmnet)
    # extract mean values
    xlab = expression(Log(lambda))
    if (sign.lambda < 0)
      xlab = paste("-", xlab, sep = "")
    plot.args = list(x = sign.lambda * log(res_obj$lambda), y = res_obj$cvm,
                     ylim = range(res_obj$cvup, res_obj$cvlo), xlab = xlab, ylab = cvobj$name,
                     type = "n")
    # draw plot
    do.call("plot", plot.args)
    # add standard error as error bars
    error_bars(sign.lambda * log(res_obj$lambda), res_obj$cvup, res_obj$cvlo,
               width = 0.01, col = "darkgrey")
    points(sign.lambda * log(res_obj$lambda), res_obj$cvm, pch = 20,
           col = "red")
    axis(side = 3, at = sign.lambda * log(res_obj$lambda), labels = paste(res_obj$nzero),
         tick = FALSE, line = -1)
    abline(v = sign.lambda * log(res_obj$lambda.min), lty = 3)
    abline(v = sign.lambda * log(res_obj$lambda.1se), lty = 3)
    invisible()
  }
}


# if using base R plot, draw error bars (copied from glmnet)
error_bars <- function (x, upper, lower, width = 0.02, ...)
{
  xlim <- range(x)
  barw <- diff(xlim) * width
  segments(x, upper, x, lower, ...)
  segments(x - barw, upper, x + barw, upper, ...)
  segments(x - barw, lower, x + barw, lower, ...)
  range(upper, lower)
}
