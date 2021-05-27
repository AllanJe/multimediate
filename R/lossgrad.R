lossgrad <- function(w,treatment,mediators,outcome,med.weights){
  
  binarymed <- is.binary(mediators)
  
  mhat.gauss <- treatment%*%w$alpha
  mhat.bin <- 1/(1+exp(-mhat.gauss))
  mhat <- mhat.gauss%*%diag((1-binarymed)) + mhat.bin%*%diag(binarymed)
  
  
  outcomehat = treatment %*% w$gamma + mediators %*% w$beta
  if (is.2(length(unique(outcome)))){
    outcomehat <- 1/(1+exp(-outcomehat))
  }
  
  wgrad <- c()
  wgrad$alpha <- t(treatment) %*% (mhat-mediators) %*% diag(med.weights) /length(outcome)
  wgrad$beta <- as.vector(t(mediators) %*% (outcomehat-outcome)) /length(outcome)
  wgrad$gamma <- as.vector(t(treatment) %*% (outcomehat-outcome))  /length(outcome)
  
  return(wgrad)
}