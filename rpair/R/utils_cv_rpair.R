#' Internal function - test if survival outcome
#'
#' Tests whether or not an outcome is a survival outcome. An outcome is assumed a survival outcome
#' if one of the following holds: (1) outcome inherits class "Surv", (2) outcome has three columns
#' (function will warn the user if column names do not match expected names), (3) outcome has two
#' columns named 'time' and 'status'. For the two-column case, if column names do not match
#' expected names, the function will assume the outcome is of type 'pairs'. An outcome with 4 or
#' more columns will cause the function to crash.
#'
#' @param y Outcome to be evaluated.
#'
#' @return Boolean value.
#'
#' @noRd
is_survival <- function(y){
  # is outcome of class Surv?
  if(inherits(y,"Surv")) return(TRUE)

  num_cols = ncol(y)

  # is outcome a vector?
  if(is.null(num_cols)) return(FALSE)

  # does outcome have three columns?
  if(num_cols == 3){
    # warn if expected column names are not present
    if(!all(colnames(y) %in% c("start", "stop", "status"))){
      warning("One or more unexpectated column names for survival class with three columns. User must ensure correct type.")
    }
    return(TRUE)
  }

  # does outcome have two columns?
  if(num_cols == 2){
    if(all(colnames(y) %in% c('time', 'status'))) return(TRUE)

    # if column names missing / don't match, assume pairs
    return(FALSE)

  }

  # y is single column
  if(num_cols==1) return(FALSE)

  # y has more than three columns
  if(num_cols > 3) stop("Outcomes with more than three columns are not supported!")


}

#' Internal function - create random folds for cross-validation
#'
#' Order non-survival type response variable and create fold splits.
#'
#' @param y Non-survival type response variable (e.g. continuous). One column only.
#' @param nfolds Number of folds to create.
#'
#' @return Data frame with columns: index, (response column), and fold.
#'
#' @noRd
get_folds <- function(y, nfolds){
  y <- cbind(index = 1:length(y), y)
  sorted <- y[order(y[,2]),]
  fold_df <- as.data.frame(create_fold_splits(sorted, nfolds))

  # handle indices with unassigned folds
  nas <- sum(is.na(fold_df$fold))
  if(nas <= nfolds){
    fold_df[is.na(fold_df$fold),"fold"] <- sample(seq(1:nfolds))[1:nas]
  }else{
    fold_df[is.na(fold_df$fold),"fold"] <- sample(seq(1:nfolds), size=nas, replace = T)
  }
  folds <- merge(sorted, fold_df, by="index")
  folds
}

#' Internal function - create random stratified folds for cross-validation
#'
#' Order survival type response variable by time then split response variable into list of event
#' groups. Create fold splits on each event group individually.
#'
#' @param S Survival type response variable. Two-columns or three-columns only.
#' @param nfolds Number of folds to create.
#'
#' @return Data frame with columns: index, (response variable columns), fold.
#'
#' @noRd
get_stratified_folds <- function(S, nfolds){
  Smat <- as.data.frame(as.matrix(S))
  num_cols <- ncol(Smat)
  Smat <- cbind(index = 1:dim(Smat)[1], Smat)

  # identify relative position of time and status columns
  if(num_cols == 2){
    time_idx <- ifelse("time" %in% colnames(Smat), which(colnames(Smat)=="time"), 2)
    status_idx <- ifelse("status" %in% colnames(Smat), which(colnames(Smat)=="status"), 3)
  }else if(num_cols == 3){
    start_idx <- ifelse("start" %in% colnames(Smat), which(colnames(Smat)=="start"), 2)
    stop_idx <- ifelse("stop" %in% colnames(Smat), which(colnames(Smat)=="stop"), 3)
    status_idx <- ifelse("status" %in% colnames(Smat), which(colnames(Smat)=="status"), 4)
    Smat <- cbind(Smat, time = Smat[,stop_idx] - Smat[,start_idx])
    time_idx = 5
  }else{
    stop("Can only create stratified folds for survival outcome with two or three columns.")
  }

  # sort survial outcome by time variable
  sorted <- Smat[order(Smat[,time_idx]),]

  # split into event groups, create fold splits for each event group
  eventlist <- split(sorted, f=sorted[,status_idx])
  fold_df_list <- lapply(eventlist, create_fold_splits, nfolds)
  fold_df <- do.call(rbind.data.frame,fold_df_list)

  # handle indices with unassigned folds
  nas <- sum(is.na(fold_df$fold))
  if(nas <= nfolds){
    fold_df[is.na(fold_df$fold),"fold"] <- sample(seq(1:nfolds))[1:nas]
  }else{
    fold_df[is.na(fold_df$fold),"fold"] <- sample(seq(1:nfolds), size=nas, replace = T)
  }
  strat_folds <- merge(sorted, fold_df, by="index")

  strat_folds
}

#' Internal function - create fold splits
#'
#' Split response variable data frame into bins with nfolds elements in each bin. Randomly assign
#' elements in each bin to a fold. If a bin contains less than nfolds elements, return NA for
#' all elements in bin.
#'
#' @param df Response variable data frame to split into fold groups.
#' @param nfolds Number of folds to create.
#'
#' @return Data frame with columns: index, fold.
#'
#' @noRd
create_fold_splits <- function(df, nfolds){
  d <- df[,"index"]
  bins <- split(d, ceiling(seq_along(d)/nfolds))
  nbins = length(bins)
  fold_list <- lapply(bins, function(x){
    if(length(x)==nfolds){
      cbind(index=x, fold=sample(seq(1,nfolds)))
    }else{
      cbind(index=x, fold=rep(NA,length(x)))
    }
  })
  folds <- do.call(rbind,fold_list)

  folds
}

#' Internal function - return loss-type string for plotting
#'
#' Return a string of the cross-validation for plotting purposes. If type.measure is deviance,
#' return the name of the loss function. If type.measure is concordance return "C-Index".
#'
#' @param type.measure Cross-validation loss function ("deviance" or "cindex").
#' @param losstype Model loss-type.
#'
#' @return Cross-validation loss string.
#'
#' @noRd
rpair_cvtype <- function(type.measure, losstype){
  if(type.measure=="deviance"){
    type.measure = losstype
    names(type.measure) = switch(losstype, lognet = "Log Loss",
                                 fishnet = "Exponential Loss",
                                 phuhnet = "Huberized Hinge Loss",
                                 psqhnet = "Squared Hinge Loss")
  }else{
    names(type.measure) <- "C-Index"
  }
  type.measure
}

#' Internal function - get min and 1se lambda values
#'
#' Find the lambda for at which the cv loss is equal to the minimum value. Also find the largest
#' value of lambda that is less than or equal to the minimum plus its corresponding standard
#' deviation.
#'
#' @param lambda Lambda sequence from master fit (excluding any values for which all folds are NA).
#' @param cvm Vector of mean loss values for each lambda.
#' @param cvsd Vector of standard deviance values for each lambda.
#' @param cvname The naem of the cv loss used. If 'cindex', reverse the direction of cvm.
#'
#' @return A list containing the lambda.min value, lambda.1se value, and the index of each.
#'
#' @noRd
getopt_cv_rpair <- function (lambda, cvm, cvsd, cvname)
{
  cvm = cvm[!is.na(lambda)]
  cvsd = cvsd[!is.na(lambda)]
  if (match(cvname, "cindex", 0)) cvm = -cvm
  cvmin = min(cvm, na.rm = TRUE)
  idmin = cvm <= cvmin
  lambda.min = max(lambda[idmin], na.rm = TRUE)
  idmin = match(lambda.min, lambda)
  semin = (cvm + cvsd)[idmin]
  id1se = cvm <= semin
  lambda.1se = max(lambda[id1se], na.rm = TRUE)
  id1se = match(lambda.1se, lambda)
  index = matrix(c(idmin, id1se), 2, 1, dimnames = list(c("min",
                                                          "1se"), "Lambda"))
  list(lambda.min = lambda.min, lambda.1se = lambda.1se, index = index)
}

#' Internal function - build prediction matrix from cross-validation folds
#'
#' Create a list of prediction matrices by appling the rpair fit object from each fold on the
#' full input data matrix, x. Predictions are calculated on the entire dataset for the sake of
#' calculating cross-validation loss using the houwelingen method.
#'
#' @param outlist List of rpair fit objects from each fold.
#' @param nlambda Number of lambda values.
#' @param lambda Lambda sequence from master fit.
#' @param x Full input data matrix.
#' @param foldid Vector of fold assignments.
#' @param alignmnet Whether the predictions in each fold are aligned according to the fraction
#'    of progress along each fold ("fraction", default) or to use the lambda values computed
#'    on the master fit to line up the data on each fold ("lambda").
#'
#' @return List of predictions matrices.
#'
#' @noRd
rpair_buildPredmat <- function(outlist, nlambda, lambda, x, foldid, alignment){
  if(alignment=="lambda") nlambda = length(lambda)

  pred_list = list()
  nfolds = max(foldid)
  nlams = double(nfolds)

  for (i in seq(nfolds)) {
    # create prediction matrix, return NA for columns in which the length of the fold
    #   lambda vector is less than nlambda
    fitobj = outlist[[i]]
    predmat = matrix(NA,nrow(x),nlambda)
    preds = switch(alignment, fraction = predict(fitobj, x),
                   lambda = predict(fitobj, x, s = lambda))
    nlami = min(ncol(preds), nlambda)
    predmat[, seq(nlami)] = preds[, seq(nlami)]

    # add rownames and column names to prediction matrix
    rn = rownames(x)
    sn = paste("s", seq(0, length = nlambda), sep = "")
    dimnames(predmat) = list(rn, sn)

    pred_list[[i]] = predmat
  }

  pred_list
}

# four rpair loss functions
exp_loss <- function(x) exp(x)
log_loss <- function(x) log(1+exp(x))
huh_loss <- function(x, delta) ((1 - -1*x - 0.5 * delta) * (-1*x <= 1 - delta) + 0.5 * (1 - -1*x)^2/delta * (-1*x <= 1) * (-1*x > 1 - delta) )
sqh_loss <- function(x) ((x >= 0) * (1+x)^2 )

#' Internal function - calculate deviance
#'
#' Uses the user-selected loss-type to calculate cross-validation loss for each lambda on each
#' fold. Can use either the standard (test pairs only) or houwelingen (all pairs - training
#' pairs) method.
#'
#' @param predmat List of prediction matrices per fold.
#' @param y Input response variable.
#' @param loss_type Loss function to use.
#' @param fids Vector of fold assignments.
#' @param delta For huberized hinge loss only, additional loss function parameter.
#' @param use_houw Whether to use the houwelingen method to calculate deviance.
#'
#' @return Data frame of nlambda x nfolds containing deviance values.
#'
#' @noRd
cv_deviance <- function(predmat, y, loss_type, fids, delta, use_houw){
  ffs <- list(exp = exp_loss, log = log_loss, huh = huh_loss, sqh = sqh_loss)
  pi_all = y_to_pairs(y)

  lo =
    sapply(seq(unique(fids)), function(i){
      if(use_houw){
        # training set samples
        j = which(fids != i)
        # all_pairs - training_set_pairs
        pir =  pi_all[!( (pi_all[,1] %in% j) & (pi_all[,2] %in% j) ),]
      }else{
        # test set samples
        j = which(fids == i)
        # test_pairs only
        pir = pi_all[( (pi_all[,1] %in% j) & (pi_all[,2] %in% j) ),]
      }

      # loss margin z i.e. loss_function(z)
      Z = predmat[[i]][pir[,1],] -  predmat[[i]][pir[,2],]
      # Z is a matrix of size:n_valid_pairs X n_lambda

      # pick loss function based on the model loss type
      ff = ffs[[loss_type]] # ffs is the list of functions

      # for each lambda calculate deviance for fold i
      if(loss_type != "phuhnet" && loss_type != "huh"){
        colMeans( ff(Z), na.rm = T )
      }else{
        colMeans( ff(Z, delta), na.rm = T )
      }
    })

  lo

}

#' Internal function - assemble concordance results
#'
#' Set up prediction matrix, response variable, and fold IDs in order to perform simultaneous
#' calculation of concordance value on all folds for each lambda value.
#'
#' @param predmat List of prediction matrices per fold.
#' @param y Input response variable.
#' @param foldid Vector of fold assignments.
#' @param use_houw Whether to use the houwelingen method to calculate concordance.
#'
#' @return Data frame with 2 rows x nlambda columns. First row is concordance, second row is
#'    variance.
#'
#' @noRd
cv_concordance <- function(predmat, y, foldid, use_houw){
  # combine predmat from each fold into single data frame
  Yh = do.call(predmat, what = rbind)
  # vector of folds indicating which fold each row of predmat belongs to
  z = unlist( lapply(seq(predmat), function(i) rep(i, nrow(predmat[[i]]))) )
  # inds of training sets for each fold
  ii = unlist(lapply(seq(predmat), function(i) foldid != i))
  # repeat response variable for each fold
  Sf = rep(y, length(predmat))

  # calculate concordance over every fold for each lambda value
  if(use_houw){
    cvfl = apply(Yh, 2, function(yh)
      tryCatch(calc_conc(Sf, yh, ii, z)["houw",], error = function(er) c(NA,NA)))
  }else{
    cvfl = apply(Yh, 2, function(yh)
      tryCatch(calc_conc(Sf, yh, ii, z)["test",], error = function(er) c(NA,NA)))
  }
  cvfl
}

#' Internal function - calculate concordance
#'
#' Calculate concordance on full data, training data, houwelingen data, and test data on lambda
#' of interest. Sample values repeated nfolds times to calculate the lambda of interest for all
#' folds at once.
#'
#' @param S Response variable repeated for each fold.
#' @param yh All sample values (repeated for each fold) for a single lambda value.
#' @param ii Boolean vector determining if index belong to training set of that fold (determined
#'    by relative position).
#' @param folds Vector of fold assignments repeated for each fold.
#'
#' @return Data frame with 4 rows (full, train, houw, and test) and 2 columns (d = concordance,
#'    v = variance).
#'
#' @noRd
calc_conc <- function(S, yh, ii, folds){

  df = data.frame(S=S, yh=yh, z=folds)
  # concordance on full data
  d0 = concordance(S~yh+strata(z), df, keepstrata = T, reverse = T)
  d0$count = rbind(d0$count)
  # concordance on train data
  d1 = concordance(S~yh+strata(z), df[ii,], keepstrata = T, reverse = T)
  d1$count = rbind(d1$count)

  #estimated concordance for houw pairs
  m = colSums(d0$count) - colSums(d1$count)
  dh = unname( (m[1]+m[3]/2)/sum(m[1:3]) )

  #--- estimate se of houw pairs
  # number of all pairs and its variance
  n0 = sum(d0$count[,1:3])
  v0 = d0$var
  d0 = d0$concordance

  # number of train pairs and variance
  n1 = sum(d1$count[,1:3])
  v1 = d1$var
  d1 = d1$concordance

  # number of houw pairs
  nh =  sum(m[1:3])
  # variance of h pairs
  vh = unname( ( n0*v0*(n1+n0) -n1*v1*n1 - n1*(d0-d1)^2 - nh*(d0-dh) )/nh )/nh

  d2 =  concordance(S~yh+strata(z), df[!ii,], reverse = T)
  cbind(
    d = c(full = 1-d0, train = 1-d1, houw = 1-dh, test = 1-d2$concordance),
    v = c(full = v0, train = v1, houw = vh, test = d2$var )
  )
}

#' Internal function - calculate statistics for each lambda
#'
#' Calculate mean and standard deviatoin of the cv loss metric across the folds for each lambda
#' value. Return only indices in which cvsd is not NA.
#'
#' @param cvfl If conc=FALSE, a data frame of nlambda x nfolds containing deviance values. Else, a
#'    data frame
#' @param lambda Lambda sequence from the master fit.
#' @param nz Vector of nonzero coefficients from the master fit.
#' @param conc Whether the cv loss type being calculated is 'concordance'. Default: FALSE.
#'
#' @return A list containing lambda, cvm, cvsd, cvup, cvlo, and nzero. See cv_rpair documentation
#'    for details.
#'
#' @noRd
cv_stats <- function (cvfl, lambda, nz, conc=F)
{
  if(conc){
    cvm = cvfl[1,]
    cvsd = sqrt(cvfl[2,])
  }else{
    cvm = rowMeans(cvfl, na.rm = T)
    cvsd = apply(cvfl, 1, sd, na.rm = T) / sqrt(ncol(cvfl))
  }
  nas = is.na(cvsd)
  if (any(nas)) {
    lambda = lambda[!nas]
    cvm = cvm[!nas]
    cvsd = cvsd[!nas]
    nz = nz[!nas]
  }
  list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvup = cvm +
         cvsd, cvlo = cvm - cvsd, nzero = nz)
}
