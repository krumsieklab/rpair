
#-------function to get comparable pairs-------#
get_riskset<- function(S){
  time <- S[,1]
  status <- S[,2]
  N = length(time)
  # for tied times
  time[status == 0] = time[status == 0]+1e-4
  dtimes <- time  
  dtimes[status == 0] = Inf
  # comparable pairs
  t(
  apply(which(outer(time, dtimes, ">"), arr.ind = T), 1, function(i){
    a = rep(0,N)
    a[i[1]] = -1
    a[i[2]] = 1
    a
  }))
}

fcvid <- function(S,k){
  if(sum(S[,2])<5) stop("no enough event!")
  
  inds1 = S[,2]==1
  inds0 = S[,2]==0
  
  k = sample(seq(k))
  ids = rep(0, length(S[,2]))
  ids[inds1] = sample(rep(k, length = sum(inds1)))
  ids[inds0] = sample(rep(k, length = sum(inds0)))
  ids
}
