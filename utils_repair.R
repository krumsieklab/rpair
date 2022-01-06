jerr<- function (fit, ...)  UseMethod("jerr")

jerr.glmnetfit <- function(fit, maxit, pmax, loss_type){
  n = fit$jerr
  
  # glmnet errors
  if (n == 0) 
    errlist = list(n = 0, fatal = FALSE, msg = "")
  else {
    outlist = {function(){
      if (n > 0) {
        if (n < 7777) 
          msg = "Memory allocation error; contact package maintainer"
        else if (n == 7777) 
          msg = "All used predictors have zero variance"
        else if (n == 10000) 
          msg = "All penalty factors are <= 0"
        else msg = "Unknown error"
        list(n = n, fatal = TRUE, msg = msg)
      }
      else if (n < 0) {
        if (n > -10000) 
          msg = paste("Convergence for ", -n, "th lambda value not reached after maxit=", maxit, 
                      " iterations; solutions for larger lambdas returned", sep = "")
        
        if (n < -10000) 
          msg = paste("Number of nonzero coefficients along the path exceeds pmax=",  pmax, " at ", -n - 10000, 
                      "th lambda value; solutions for larger lambdas returned", sep = "")
        
        list(n = n, fatal = FALSE, msg = msg)
      }
    }}()
    errlist = switch(loss_type, 
                     log = {function(){
                       if (n < -20000) 
                         outlist$msg = paste("Max(p(1-p),1.0e-6 at ", -n - 20000, "th value of lambda; solutions for larger values of lambda returned")
                       if (outlist$msg != "Unknown error") 
                         return(outlist)
                       if ( ((8000 < n) & (n < 9000)) | ((9000 < n) & (n < 10000)) )
                         msg = paste("Null probability for concordance < 1.0e-5")
                       else msg = "Unknown error"
                       list(n = n, fatal = TRUE, msg = msg)
                     }}(),  
                     exp = {function(){
                       if (outlist$msg != "Unknown error") 
                         return(outlist)
                       if (n == 8888)  msg = "NA error"
                       else if (n == 9999) msg = "No positive observation weights"
                       else msg = "Unknown error"
                       list(n = n, fatal = TRUE, msg = msg)
                     }}())
    
    names(errlist) = c("n", "fatal", "msg")
    errlist$msg = paste("from glmnet Fortran code (error code ",  n, "); ", errlist$msg, sep = "")
  }
  return(errlist)
}

jerr.gcdnetfit <- function(fit, maxit, pmax, loss_type){

  n = fit$jerr
  # gcdnet errors
  if (n == 0) 
    msg <- ""
  if (n > 0) {
    if (n < 7777) 
      msg <- "Memory allocation error"
    if (n == 7777) 
      msg <- "All used predictors have zero variance"
    if (n == 10000) 
      msg <- "All penalty factors are <= 0"
    n <- 1
    msg <- paste("in gcdnet fortran code -", msg)
  }
  if (n < 0) {
    if (n > -10000) 
      msg <- paste("Convergence for ", -n, "th lambda value not reached after maxit=", 
                   maxit, " iterations; solutions for larger lambdas returned", 
                   sep = "")
    if (n < -10000) 
      msg <- paste("Number of nonzero coefficients along the path exceeds pmax=", 
                   pmax, " at ", -n - 10000, "th lambda value; solutions for larger lambdas returned", 
                   sep = "")
    n <- -1
    msg <- paste("from gcdnet fortran code -", msg)
  }
  
  list(n = n, msg = msg)
  
}

fix_lam<- function (lam) {
  if (length(lam) > 2) {
    llam = log(lam)
    lam[1] = exp(2 * llam[2] - llam[3])
  }
  lam
}

getcoef<- function (fit, ...)  UseMethod("getcoef")

getcoef.glmnetfit <- function (fit, nvars, pmax, vnames) {
  
  lmu = fit$lmu
  if (lmu < 1) {
    warning("an empty model has been returned; probably a convergence issue")
    coefob = list(beta = zeromat(nvars, as.integer(1), vnames, "s0"), df = 0, dim = c(nvars, 1), lambda = Inf)
    return(coefob)
  }
  nin = fit$nin[seq(lmu)]
  ninmax = max(nin)
  lam = fit$alm[seq(lmu)]
  stepnames = paste("s", seq(lmu) - 1, sep = "")
  dd = c(nvars, lmu)
  if (ninmax > 0) {
    ca = matrix(fit$ca[seq(pmax * lmu)], pmax, lmu)[seq(ninmax), , drop = FALSE]
    df = apply(abs(ca) > 0, 2, sum)
    ja = fit$ia[seq(ninmax)]
    oja = order(ja)
    ja = rep(ja[oja], lmu)
    ia = cumsum(c(1, rep(ninmax, lmu)))
    beta = Matrix::drop0( new("dgCMatrix", Dim = dd, Dimnames = list(vnames, stepnames), 
                              x = as.vector(ca[oja, ]), p = as.integer(ia - 1), i = as.integer(ja - 1)) )
  }
  else {
    beta = zeromat(nvars, lmu, vnames, stepnames)
    df = rep(0, lmu)
  }
  
  list(beta = beta, df = df, dim = dd, lambda = lam)
}

getcoef.gcdnetfit <- function (fit, nvars, pmax, vnames) { #maxit, 
  
  nalam <- fit$nalam
  nbeta <- fit$nbeta[seq(nalam)]
  nbetamax <- max(nbeta)
  lam <- fit$alam[seq(nalam)]
  stepnames <- paste("s", seq(nalam) - 1, sep = "")
  
  # errmsg <- err(fit$jerr, maxit, pmax)
  # 
  # switch(paste(errmsg$n), 
  #        `1` = stop(errmsg$msg, call. = FALSE), 
  #        `-1` = print(errmsg$msg, call. = FALSE))
  # 
  
  
  dd <- c(nvars, nalam)
  if (nbetamax > 0) {
    beta <- matrix(fit$beta[seq(pmax * nalam)], pmax, nalam)[seq(nbetamax), , drop = FALSE]
    
    df <- apply(abs(beta) > 0, 2, sum)
    ja <- fit$ibeta[seq(nbetamax)]
    oja <- order(ja)
    ja <- rep(ja[oja], nalam)
    ibeta <- cumsum(c(1, rep(nbetamax, nalam)))
    beta <- new("dgCMatrix", 
                Dim = dd, 
                Dimnames = list(vnames, stepnames), 
                x = as.vector(beta[oja, ]), 
                p = as.integer(ibeta - 1), 
                i = as.integer(ja - 1))
  }
  else {
    beta <- zeromat(nvars, nalam, vnames, stepnames)
    df <- rep(0, nalam)
  }
  
  list( beta = beta, df = df, dim = dd, lambda = lam)
}

zeromat <- function (nvars, nalam, vnames, stepnames) {
  ca <- rep(0, nalam)
  ia <- seq(nalam + 1)
  ja <- rep(1, nalam)
  dd <- c(nvars, nalam)
  new("dgCMatrix", Dim = dd, Dimnames = list(vnames, stepnames), 
      x = as.vector(ca), p = as.integer(ia - 1), i = as.integer(ja - 1))
}

getmin <- function (lambda, cvm, cvsd) {
  cvmin = min(cvm, na.rm = TRUE)
  idmin = cvm <= cvmin
  lambda.min = max(lambda[idmin], na.rm = TRUE)
  idmin = match(lambda.min, lambda)
  semin = (cvm + cvsd)[idmin]
  idmin = cvm <= semin
  lambda.1se = max(lambda[idmin], na.rm = TRUE)
  list(lambda.min = lambda.min, lambda.1se = lambda.1se)
}

