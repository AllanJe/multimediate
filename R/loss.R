loss <- function(w,treatment,mediators,outcome,med.weights){
  
  K=dim(mediators)[2]
  P=dim(treatment) [2]
  binarymed <- is.binary(mediators)
  mhat <- treatment %*% w$alpha
  zhat = treatment %*% w$gamma + mediators %*% w$beta
  
  
  loss<- sum(mapply(function(x,y,bin) variableLoss(x,y,bin),mhat,mediators,binarymed)*med.weights) + 
    variableLoss(zhat,outcome,is.2(length(unique(outcome))))
  
  return(loss/length(outcome))
}