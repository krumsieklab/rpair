#' ADD TITLE
#'
#' ADD DESCRIPTION
#'
#' @param x
#' @param y
#' @param loss_type Loss function to use. One of c("log", "exp").
#' @param alpha Default: 1.0 (lasso).
#' @param nlambda Default: 100.
#' @param lambda.min.ratio Default: ifelse(nobs0 <= nvars, 1e-2, 1e-4)
#' @param lambda Default: NULL.
#' @param thresh
#' @param
#'
#' @author mubu
#'
#' @export
repair_gloss<-
  function( x, # data
            y, # outcome
            loss_type=c("log","exp"), # equivalent to 'family'
            weights, # x
            offset=NULL, # x
            alpha=1.0, # yes
            nlambda=100, # yes
            lambda.min.ratio=ifelse(nobs0 <= nvars,1e-2,1e-4), # yes
            lambda=NULL, # yes
            standardize=FALSE, # maybe
            thresh=1e-7, # yes
            dfmax=nvars+1, #yes
            pmax=min(dfmax*2+20,nvars), # yes
            exclude, # x
            penalty.factor=rep(1,nvars), # yes
            lower.limits=-Inf, # ?? GGM
            upper.limits=Inf, # ??
            maxit=100000, # yes
            type.logistic=c("Newton","modified.Newton")#, # yes
            #trace.it = 0 # yes
            ){
    print(lambda[1:5])
    # print(lambda.min.ratio)
    nobs0 = nrow(x)
    # if( is.matrix(y) & !(identical(colnames(y),c("time","status")) | identical(colnames(y),c("start", "stop", "status")) ) ){
    #   if( any(!(y[,1] %in% seq(nobs0))) | any(!(y[,2] %in% seq(nobs0))) )
    #     stop("some observations in the pairs don't exsist")
    # }

    cp = y_to_pairs(y)
    ncp = nrow(cp)

    ### Prepare all the generic arguments, then hand off to loss_type functions
    loss_type=match.arg(loss_type)
    if(alpha>1){
      warning("alpha >1; set to 1")
      alpha=1
    }
    if(alpha<0){
      warning("alpha<0; set to 0")
      alpha=0
    }
    alpha=as.double(alpha)

    this.call=match.call() # KC: What is this for?
    nlam=as.integer(nlambda) # KC: Why rename this?
    np=c(ncp,ncol(x))

    ###check dims
    if(is.null(np)|(np[2]<=1))stop("x should be a matrix with 2 or more columns")
    nobs=as.integer(np[1])
    if(missing(weights))weights=rep(1,nobs)
    else if(length(weights)!=nobs)stop(paste("number of elements in weights (",length(weights),") not equal to the number of rows of cp (",nobs,")",sep=""))
    nvars=as.integer(np[2])

    vnames=colnames(x)
    if(is.null(vnames))vnames=paste("V",seq(nvars),sep="")
    ne=as.integer(dfmax)
    nx=as.integer(pmax)
    if(missing(exclude))exclude=integer(0)
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

    # -- here will need special treatment later
    ###check on limits
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
      #message("for limits==0, needs glmnet::glmnet.control()")
      require(glmnet)
      ###Bounds of zero can mess with the lambda sequence and fdev; ie nothing happens and if fdev is not
      ###zero, the path can stop
      fdev=glmnet::glmnet.control()$fdev
      if(fdev!=0) {
        glmnet::glmnet.control(fdev=0)
        on.exit(glmnet::glmnet.control(fdev=fdev))
      }
    }
    storage.mode(cl)="double"
    ### end check on limits

    isd=as.integer(standardize)
    thresh=as.double(thresh)

    cat("\nacaba: ")
    print(is.null(lambda))

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

    fit <-
      switch(loss_type,
             log = plognet(x,cp,weights,offset,alpha,nobs,nvars,jd,vp,cl,ne,nx,nlam,flmin,ulam,thresh,isd,vnames,maxit,kopt),
             exp = pfishnet(x,cp,weights,offset,alpha,nobs,nvars,jd,vp,cl,ne,nx,nlam,flmin,ulam,thresh,isd,vnames,maxit))


    if(is.null(lambda))fit$lambda=fix_lam(fit$lambda) #glmnet:::fix.lam(fit$lambda)##first lambda is infinity; changed to entry point
    fit$call=this.call
    fit$loss=loss_type
    fit$nobs=nobs

    class(fit)=c(class(fit),"glmnet")
    fit
  }

plognet <-
  function( x,cp,weights,offset,
            alpha,nobs,nvars,jd,vp,cl,ne,nx,nlam,flmin,ulam,thresh,isd,
            vnames,maxit,kopt, rk = nrow(cp)%/%2
            ){
  print("HELLO")
  # modify cp accordingly
  cp = rbind( cp[1:rk,2:1], cp[-(1:rk),])
  y = c(rep(1,rk),rep(0, nrow(cp)-rk))
  y = cbind(c0=1-y, c1 =y)
  # ---

  maxit=as.integer(maxit)
  nc=2
  #y=diag(nc)[as.numeric(y),]
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

  #loss_type="log"
  nc=as.integer(1) # for calling lognet
  y=y[,c(2,1)]#fortran lognet models the first column; we prefer the second (for 0/1 data)

  storage.mode(y)="double"
  if(is.null(offset)){
    offset=y*0 #keeps the shape of y
    is.offset=FALSE}
  else{
    offset=as.matrix(offset)
    if(!is.null(o))offset=offset[o,,drop=FALSE]# we might have zero weights
    do=dim(offset)
    if(do[[1]]!=nobs)stop("offset should have the same number of values as pairs",call.=FALSE)
    if((do[[2]]==1)&(nc==1))offset=cbind(offset,-offset)
    storage.mode(offset)="double"
    is.offset=TRUE
  }

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
    errmsg=jerr(fit,maxit,pmax=nx,"log") #glmnet::jerr(fit$jerr,maxit,pmax=nx,"binomial")
    if(errmsg$fatal)stop(errmsg$msg,call.=FALSE)
    else warning(errmsg$msg,call.=FALSE)
  }

  outlist=getcoef(fit,nvars,nx,vnames)

  dev=fit$dev[seq(fit$lmu)]
  outlist=c(outlist,list(npasses=fit$nlp, jerr=fit$jerr,dev.ratio=dev,nulldev=fit$nulldev,offset=is.offset))

  # class(outlist)="lognet"
  class(outlist)="repair"
  outlist
}

pfishnet <-
  function( x,cp, weights,offset,
            alpha,nobs,nvars,jd,vp,cl,ne,nx,nlam,flmin,ulam,thresh,isd,
            vnames,maxit
            ){
  print("hi from fishnet")

  maxit=as.integer(maxit)
  weights=as.double(weights)

  if(is.null(offset)){
    offset=double(nrow(cp)) #keeps the shape of y
    is.offset=FALSE}
  else{
    storage.mode(offset)="double"
    is.offset=TRUE
  }
  #loss_type="exp"
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

  outlist=getcoef(fit,nvars,nx,vnames) #glmnet::getcoef(fit,nvars,nx,vnames)
  dev=fit$dev[seq(fit$lmu)]
  outlist=c(outlist,list(npasses=fit$nlp,jerr=fit$jerr,dev.ratio=dev,nulldev=fit$nulldev,offset=is.offset))
  # class(outlist)="fishnet"

  class(outlist)="repair"
  outlist
}

