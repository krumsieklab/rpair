#' Fit GLM for pairwise survival ranking models
#'
#' Fit a generalized linear model for pairwise survival ranking models using either logistic loss (concordance
#' regression) or exponential loss. Uses a modified implementation of \link[glmnet]{glmnet}. Refer to glmnet
#' documentation for further details.
#'
#' @param x Input matrix of dimension nobs x nvars; each row is an observation vector.
#' @param y Pairwise ranking analysis response variable. This function supports three types of input types:
#'    (1) continuous values, (2) survival data, and (3) ranked pairs.
#' @param loss_type Loss function to use. Equivalent to "family" in glmnet package. One of c("log", "exp").
#' @param alpha The elasticnet mixing parameter. See \link[glmnet]{glmnet} for details. Default: 1.0 (lasso penalty).
#' @param nlambda Number of lambda values to evaluate. Default: 100.
#' @param lambda.min.ratio Smallest value for lambda as a fraction of lambda.max, the (data derived) entry value.
#'   See \link[glmnet]{glmnet} for details. Default: ifelse(nobs0 <= nvars, 1e-2, 1e-4)
#' @param lambda A user supplied lambda sequence. Overrides the typical usage in which a lambda sequence is computed
#'    using nlambda and lambda.min.ratio. Provide a decreasing sequence of lambda values with at least 2 entries.
#'    Default: NULL.
#' @param standardize Logical flag for x variable standardization. Default: FALSE.
#' @param thresh Convergence threshold for coordinate descent. Default: 1e-7.
#' @param dfmax Limit the maximum number of variables in the model. Default: nvars+1.
#' @param pmax Limit the maximum number of variables that can be nonzero. Default: min(dfmax*2+20, nvars).
#' @param penalty.factor Vector of penalty factors to apply to each coefficient. Default: rep(1, nvars).
#' @param maxit Maximum number of passes over the data for all lambda values. Default: 100000.
#' @param type.logistic Only used by logistic loss. One of c("Newton", "modified.Newton"). If "Newton" then the
#'    exact hessian is used, while "modified.Newton" uses an upper-bound on the hessian and can be faster.
#'    Default: "Newton".
#' @param use_glmnet Whether to use the glmnet function glmnet.control for handling changes needed for upper
#'    and lower limits. Only use if glment is installed. Uses rpair_control if FALSE. Default: FALSE.
#'
#' @return An object with S3 class \code{"rpair"}, "*", where "*" is \code{"lognet"} or \code{"fishnet"}. Contains
#' the following attributes:
#'  \itemize{
#'    \item{beta - a nvars x length(lambda) matrix of coefficients, stored in sparse column format}
#'    \item{df - the number of nonzero coefficients for each value of lambda}
#'    \item{dim - dimension of coefficient matrix}
#'    \item{lambda - The actual sequence of lambda values used}
#'    \item{npasses - total passes over the data summed over all lambda values}
#'    \item{jerr - error flag, for warnings and errors (largely for internal debugging)}
#'    \item{dev.ratio - The fraction of (null) deviance explained}
#'    \item{nulldev - the null deviance (per observation)}
#'    \item{call - the call that produced the object}
#'    \item{loss - the loss function used}
#'    \item{nobs - the number of observations}
#'  }
#'
#'
#' @examples
#' efit = rpair_gloss(surv_x, surv_cp, pmax = 50, loss_type = "exp")
#' lfit = rpair_gloss(surv_x, surv_cp, pmax = 50, loss_type = "log")
#'
#' @author mubu, KC
#'
#' @export
rpair_gloss<-
  function( x,
            y,
            loss_type=c("exp","log"),
            alpha=1.0,
            nlambda=100,
            lambda.min.ratio=ifelse(nobs0 <= nvars,1e-2,1e-4),
            lambda=NULL,
            standardize=FALSE,
            thresh=1e-7,
            dfmax=nvars+1,
            pmax=min(dfmax*2+20,nvars),
            penalty.factor=rep(1,nvars),
            maxit=100000,
            type.logistic=c("Newton","modified.Newton"),
            use_glmnet = FALSE
  ){

    # return call with fitted object
    this.call=match.call()

    # previously user-defined variables - for now not supporting the ability for user to set them
    exclude=integer(0)
    lower.limits=-Inf
    upper.limits=Inf
    # offset - defined in loss functions
    # weights - defined with generic parameters

    ### Prepare input and outcome variables
    # for evaluating lambda.min.ratio only
    nobs0 = nrow(x)

    # generate comparable pairs
    cp = y_to_pairs(y)
    ncp = nrow(cp)

    #check matrix dimensions
    np=c(ncp,ncol(x))
    if(is.null(np)|(np[2]<=1))stop("x should be a matrix with 2 or more columns")

    nobs=as.integer(np[1])
    nvars=as.integer(np[2])

    vnames=colnames(x)
    if(is.null(vnames))vnames=paste("V",seq(nvars),sep="")

    ### Prepare all the generic arguments, then hand off to loss_type functions
    ## unmodified parameters
    weights=rep(1,nobs)
    loss_type=match.arg(loss_type)
    ne=as.integer(dfmax)
    nx=as.integer(pmax)
    isd=as.integer(standardize)
    thresh=as.double(thresh)

    ## parameters with checks / conditions
    if(alpha>1){
      warning("alpha >1; set to 1")
      alpha=1
    }
    if(alpha<0){
      warning("alpha<0; set to 0")
      alpha=0
    }
    alpha=as.double(alpha)

    if(any(penalty.factor==Inf)){
      exclude=c(exclude,seq(nvars)[penalty.factor==Inf])
      exclude=sort(unique(exclude))
    }
    if(length(exclude)>0){
      jd=match(exclude,seq(nvars),0)
      if(!all(jd>0))stop("Some excluded variables out of range")
      penalty.factor[jd]=1 #ow can change lambda sequence
      jd=as.integer(c(length(jd),jd))
    }else jd=as.integer(0)
    vp=as.double(penalty.factor)

    nlam=as.integer(nlambda)
    if(is.null(lambda)){
      if(lambda.min.ratio>=1)stop("lambda.min.ratio should be less than 1")
      flmin=as.double(lambda.min.ratio)
      ulam=double(1)
    }
    else{
      flmin=as.double(1)
      if(any(lambda<0))stop("lambdas should be non-negative")
      ulam=as.double(rev(sort(lambda)))
      nlam=as.integer(length(lambda))
    }

    kopt=switch(match.arg(type.logistic),
                "Newton"=0,#This means to use the exact Hessian
                "modified.Newton"=1 # Use the upper bound
    )
    kopt=as.integer(kopt)

    # -- here will need special treatment later
    ### Check lower and upper limits
    if(any(lower.limits>0)){stop("Lower limits should be non-positive")}
    if(any(upper.limits<0)){stop("Upper limits should be non-negative")}
    lower.limits[lower.limits==-Inf]= -9.9e+35
    upper.limits[upper.limits==Inf]= 9.9e+35
    if(length(lower.limits)<nvars){
      if(length(lower.limits)==1)lower.limits=rep(lower.limits,nvars)else stop("Require length 1 or nvars lower.limits")
    }
    else lower.limits=lower.limits[seq(nvars)]
    if(length(upper.limits)<nvars){
      if(length(upper.limits)==1)upper.limits=rep(upper.limits,nvars)else stop("Require length 1 or nvars upper.limits")
    }
    else upper.limits=upper.limits[seq(nvars)]
    cl=rbind(lower.limits,upper.limits)
    # where glmnet package is needed
    if(any(cl==0)){
      # Bounds of zero can mess with the lambda sequence and fdev; ie nothing happens and if fdev is not
      #    zero, the path can stop
      if(use_glmnet){
        require(glmnet)
        fdev=glmnet::glmnet.control()$fdev
        if(fdev!=0) {
          glmnet::glmnet.control(fdev=0)
          on.exit(glmnet::glmnet.control(fdev=fdev))
        }
      }else{
        fdev=rpair:::rpair_control_test()$fdev
        if(fdev!=0){
          rpair:::rpair_control_test(fdev=0)
          on.exit(rpair:::rpair_control_test(fdev=fdev))
        }
      }
    }
    storage.mode(cl)="double"
    ### end check on limits
    # glmnet.control function injects parameters into fortran code
    # fix this at the end if there is time

    ### Fit model
    fit <-
      switch(loss_type,
             log = plognetfit(x,cp,weights,alpha,nobs,nvars,jd,vp,cl,ne,nx,nlam,flmin,ulam,thresh,isd,vnames,maxit,kopt),
             exp = pfishnetfit(x,cp,weights,alpha,nobs,nvars,jd,vp,cl,ne,nx,nlam,flmin,ulam,thresh,isd, vnames,maxit))


    if(is.null(lambda))fit$lambda=fix_lam(fit$lambda)
    fit$call=this.call
    fit$loss=loss_type
    fit$nobs=nobs

    class(fit)=c(class(fit),"rpair")
    fit
  }


#' Logistic Loss (Concordance Regression) Function
#'
#' Helper function to call the Fortran implemented logistic loss (concordance regression) algorithm.
#'
#' @param x Input matrix of dimension nobs x nvars; each row is an observation vector.
#' @param cp Survival data outcome as comparable pairs.
#' @param weights Observation weights. Currently only supports 1 for each observation.
#' @param alpha Elasticnet mixing parameter (between 0-1).
#' @param nobs Number of samples (first dimension of x matrix).
#' @param nvars Number of features (second dimension of x matrix).
#' @param jd A vector of features to be excluded.
#' @param vp Vector of penalty factors to apply to each coefficient.
#' @param cl A matrix of lower and upper limit values for each coefficient.
#' @param ne Limit the maximum number of variables in the model.
#' @param nx Limit the maximum number of variables that can be nonzero.
#' @param nlam Number of lambda values to evaluate.
#' @param flmin If lambda sequence not provided, this is the lambda.min.ratio. If it is, then it is equal to 1.
#' @param ulam If provided, this is the user lambda sequence. If not, this is 1.
#' @param thresh Convergence threshold for coordinate descent.
#' @param isd Integer value for rpair_gloss 'standardize' argument.
#' @param vnames Column names of the input matrix.
#' @param maxit Maximum number of passes over the data for all lambda values.
#' @param kopt Integer value of type.logistic - 0 for Newton, 1 for modified.Newton.
#' @param rk Half the number of pairs.
#'
#' @author mubu, KC
#'
#' @noRd
plognetfit <-
  function( x,cp,weights,
            alpha,nobs,nvars,jd,vp,cl,ne,nx,nlam,flmin,ulam,thresh,isd,
            vnames,maxit,kopt, rk = nrow(cp)%/%2
  ){
    # modify cp accordingly
    cp = rbind( cp[1:rk,2:1], cp[-(1:rk),])
    y = c(rep(1,rk),rep(0, nrow(cp)-rk))
    y = cbind(c0=1-y, c1 =y)
    # ---

    maxit=as.integer(maxit)
    nc=2
    #Check for size limitations
    maxvars=.Machine$integer.max/(nlam*nc)
    if(nx>maxvars)
      stop(paste("Integer overflow; 2*num_lambda*pmax should not exceed .Machine$integer.max. Reduce pmax to be below", trunc(maxvars)))

    if(!missing(weights)) y=y*weights
    ### check if any rows of y sum to zero, and if so deal with them
    weights=drop(y%*%rep(1,nc))
    o=weights>0
    if(!all(o)){ #subset the data
      y=y[o,]
      x=x[o,,drop=FALSE]
      nobs=sum(o)
    }else o=NULL

    nc=as.integer(1) # for calling lognet
    y=y[,c(2,1)] #fortran lognet models the first column; we prefer the second (for 0/1 data)

    storage.mode(y)="double"
    offset=y*0 #keeps the shape of y

    fit <- .Fortran( "plognet",
                     parm=alpha,nrow(x),nvars,nc,as.double(x),y,offset,jd,vp,cl,ne,nx,nlam,flmin,ulam,thresh,isd,maxit,kopt,
                     lmu=integer(1),
                     a0=double(nlam*nc),
                     ca=double(nx*nlam*nc),
                     ia=integer(nx),
                     nin=integer(nlam),
                     nulldev=double(1),
                     dev=double(nlam),
                     alm=double(nlam),
                     nlp=integer(1),
                     nobs, as.integer(cp),
                     jerr=integer(1)
    )
    class(fit) = c("glmnetfit",class(fit))

    if(fit$jerr!=0){
      errmsg=jerr(fit,maxit,pmax=nx,"log")
      if(errmsg$fatal)stop(errmsg$msg,call.=FALSE)
      else warning(errmsg$msg,call.=FALSE)
    }

    outlist=getcoef(fit,nvars,nx,vnames)

    dev=fit$dev[seq(fit$lmu)]
    # KC: add 'type' to this list
    outlist=c(outlist,list(npasses=fit$nlp, jerr=fit$jerr,dev.ratio=dev,nulldev=fit$nulldev))

    class(outlist)="lognet"
    outlist
  }


#' Exponential Loss Function
#'
#' Helper function to call the Fortran implemented exponential loss algorithm.
#'
#' @param x Input matrix of dimension nobs x nvars; each row is an observation vector.
#' @param cp Comparable pairs of the form: two column matrix in which numbers correspond to sample indices;
#'   each row is a pair in which the sample in the left hand column has a greater survival time than the sample in
#'   the right hand column.
#' @param weights Observation weights. Currently only supports 1 for each observation.
#' @param alpha Elasticnet mixing parameter (between 0-1).
#' @param nobs Number of samples (first dimension of x matrix).
#' @param nvars Number of features (second dimension of x matrix).
#' @param jd A vector of features to be excluded.
#' @param vp Vector of penalty factors to apply to each coefficient.
#' @param cl A matrix of lower and upper limit values for each coefficient.
#' @param ne Limit the maximum number of variables in the model.
#' @param nx Limit the maximum number of variables that can be nonzero.
#' @param nlam Number of lambda values to evaluate.
#' @param flmin If lambda sequence not provided, this is the lambda.min.ratio. If it is, then it is equal to 1.
#' @param ulam If provided, this is the user lambda sequence. If not, this is 1.
#' @param thresh Convergence threshold for coordinate descent.
#' @param isd Integer value for rpair_gloss 'standardize' argument.
#' @param vnames Column names of the input matrix.
#' @param maxit Maximum number of passes over the data for all lambda values.
#'
#' @author mubu, KC
#'
#' @noRd
pfishnetfit <-
  function( x,cp, weights,
            alpha,nobs,nvars,jd,vp,cl,ne,nx,nlam,flmin,ulam,thresh,isd,
            vnames,maxit
  ){

    maxit=as.integer(maxit)
    weights=as.double(weights)
    offset=double(nrow(cp)) #keeps the shape of y

    fit <- .Fortran( "pfishnet",
                     parm=alpha,nrow(x),nvars,as.double(x),#y,
                     offset,weights,jd,vp,cl,ne,nx,nlam,flmin,ulam,thresh,isd,#intr,
                     maxit,
                     lmu=integer(1),
                     a0=double(nlam),
                     ca=double(nx*nlam),
                     ia=integer(nx),
                     nin=integer(nlam),
                     nulldev=double(1),
                     dev=double(nlam),
                     alm=double(nlam),
                     nlp=integer(1),
                     nobs, as.integer(cp),
                     jerr=integer(1)
    )
    class(fit) = c("glmnetfit",class(fit))

    if(fit$jerr!=0){
      errmsg=jerr(fit,maxit,pmax=nx,"exp")
      if(errmsg$fatal)stop(errmsg$msg,call.=FALSE)
      else warning(errmsg$msg,call.=FALSE)
    }

    outlist=getcoef(fit,nvars,nx,vnames)
    dev=fit$dev[seq(fit$lmu)]
    outlist=c(outlist,list(npasses=fit$nlp,jerr=fit$jerr,dev.ratio=dev,nulldev=fit$nulldev))#,offset=is.offset))

    class(outlist)="fishnet"
    outlist
  }

