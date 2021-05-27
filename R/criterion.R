criterion <- function(w,treatment,mediators,outcome,mu,groups,tau){
  K=dim(mediators)[2]
  omega <- 0
  for (k in 1:K){
    omega <- omega + norm(as.matrix( c(w$alpha[groups$alpha==k],w$beta[groups$beta==k]) ),type="F")
  }
  return(f(w,treatment,mediators,outcome,tau)+mu*omega)
}
