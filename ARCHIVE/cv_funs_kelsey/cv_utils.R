# predict.cv_rpair
predict.cv_rpair <- function (object, newx, s = c("lambda.1se", "lambda.min"), ...) 
{
  if (is.numeric(s)) 
    lambda = s
  else if (is.character(s)) {
    s = match.arg(s)
    lambda = object[[s]]
    names(lambda) = s
  }
  else stop("Invalid form for s")
  predict(object$rpair.fit, newx, s = lambda, ...)
}

# coef.cv_rpair
coef.cv_rpair <- function (object, s = c("lambda.1se", "lambda.min"), ...) 
{
  if (is.numeric(s)) 
    lambda = s
  else if (is.character(s)) {
    s = match.arg(s)
    lambda = object[[s]]
  }
  else stop("Invalid form for s")
  coef(object$rpair.fit, s = lambda, ...)
}

# plot.cv_rpair
# rearrange geom_points so dots over errorbars
# add text to top of plot with nzero values
# replace "poisson deviance" with "log loss"
# check glmnet implementation to see how se is calculated
plot.cv_rpair <- function (cvobj, sign.lambda = 1, ggplot=F) #, ...) 
{
  if(ggplot){
    plot_df <- data.frame(cvm = cvobj$cvm, cvup = cvobj$cvup, cvlo = cvobj$cvlo, nzero = cvobj$nzero) %>% 
      mutate(log_lambda = log(cvobj$lambda)) %>%
      filter(!is.na(log_lambda))
    
    plot_df %>% ggplot(aes(x=log_lambda,y=cvm)) + 
      geom_point(color="red") +
      xlab("Log(lambda)") +
      ylab(cvobj$name) +
      theme_minimal() +
      geom_errorbar(aes(ymin=cvlo, ymax=cvup)) +
      geom_vline(xintercept = log(cvobj$lambda.min), linetype="dashed") +
      geom_vline(xintercept = log(cvobj$lambda.1se), linetype="dashed")
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
