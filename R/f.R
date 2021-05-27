f <- function(w,treatment,mediators,outcome,tau){

  #print(dim(treatment))
  #print(dim(w$alpha))
  K=dim(mediators)[2]
  P=dim(treatment) [2]
  nuhat = treatment %*% w$alpha
  zhat = treatment %*% w$gamma + mediators %*% w$beta

  f=0
  for (k in 1:K){
    if (is.binary(mediators[,k])){
      f= f-sum(mediators[,k]*nuhat[,k])+sum(log(1+exp(nuhat[,k])))*tau
    }
    else{
      f=f+ .5*sqrt(sum((nuhat[,k]-mediators[,k])^2))*tau
    }
  }
  if (is.binary(outcome)){
    f= f-sum(outcome*zhat)+sum(log(1+exp(zhat)))
  }
  else{
    f=f+ .5*sqrt(sum((zhat-outcome)^2))
  }

  return(f/length(outcome))
}
