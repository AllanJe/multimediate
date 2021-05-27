proximalapprox <- function(u,mu,groups){

  K=length(u$beta)

  for(k in 1:K){
    knorm=norm(as.matrix( c(u$alpha[groups$alpha==k],u$beta[groups$beta==k]) ),type="F")
    coeff = ifelse(1-mu/knorm>0, 1-mu/knorm, 0)   #coefficient par lequel multiplier les elements du groupe k
    u$alpha[groups$alpha==k] <- u$alpha[groups$alpha==k]*coeff
    u$beta[groups$beta==k] <- u$beta[groups$beta==k]*coeff
  }
  return(u)
}
