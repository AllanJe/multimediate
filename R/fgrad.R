fgrad <- function(w,treatment,mediators,outcome,tau){

  mhat <- treatment%*%w$alpha
  for (k in 1:dim(mediators)[2]){
    if (is.binary(mediators[,k])){
      mhat[,k] <- 1/(1+exp(-mhat[,k]))
    }
  }

  outcomehat = treatment %*% w$gamma + mediators %*% w$beta
  if (is.binary(outcome)){
    outcomehat <- 1/(1+exp(-outcomehat))
  }

  wgrad <- c()
  wgrad$alpha <- t(treatment) %*% (mhat-mediators) /length(outcome) *tau
  wgrad$beta <- as.vector(t(mediators) %*% (outcomehat-outcome)) /length(outcome)
  wgrad$gamma <- as.vector(t(treatment) %*% (outcomehat-outcome))  /length(outcome)

  return(wgrad)
}
