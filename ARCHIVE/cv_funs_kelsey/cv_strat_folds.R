# cross-validation stratified folds
get_stratified_folds <- function(S, nfolds){
  Smat <- as.data.frame(as.matrix(S))
  Smat <- cbind(index = 1:dim(Smat)[1], Smat)
  sorted <- Smat[order(Smat$time),]
  eventlist <- split(sorted, f=sorted$status)
  fold_df_list <- lapply(eventlist, create_fold_splits, nfolds)
  fold_df <- do.call(rbind.data.frame,fold_df_list)
  fold_df[is.na(fold_df$fold),"fold"] <- sample(seq(1:nfolds))[1:sum(is.na(fold_df$fold))]
  strat_folds <- merge(sorted, fold_df, by="index")
  strat_folds
}

# avoid dplyr dependency
# do without calling extra package
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


ks = 100

# number of features
n1 = 20 # n informative features
n2 = 40 # n noisy features

set.seed(42)

# coefficients
b = c( exp(runif(n1)), seq(n2)*0) %>% round(2)

x = matrix( rnorm((n1+n2)*100), nrow = ks ) %>% scale
# survival outcome
S = Surv(exp( scale(x %*% b)), sample(c(F,T),ks, T) )
colnames(S) = c("time", "status")

fids <- get_stratified_folds(S, nfolds=5)
table(fids$fold)

fids <- get_stratified_folds(S, nfolds=4)
table(fids$fold)

fids <- get_stratified_folds(S, nfolds=6)
table(fids$fold)
