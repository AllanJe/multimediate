oneIteration <- function(w,treatment,mediators,outcome,groups,lambda,L0,eta,tau=1){

  L <- L0

  u <- c()
  wgrad <- fgrad(w,treatment=treatment,mediators=mediators,outcome=outcome,tau=tau)
  u$alpha <- w$alpha - 1/L*wgrad$alpha
  u$beta <- w$beta - 1/L*wgrad$beta
  u$gamma <- w$gamma - 1/L*wgrad$gamma

  wL <- proximalapprox(u,mu=lambda/L,groups=groups)

  test <-  f(w,treatment,mediators,outcome,tau=tau) - f(wL,treatment,mediators,outcome,tau=tau) +
    sum(wgrad$alpha*(wL$alpha-w$alpha))+
    sum(wgrad$beta*(wL$beta-w$beta))+
    sum(wgrad$gamma*(wL$gamma-w$gamma))+
    L/2*norm((wL$alpha-w$alpha),type="F")^2 +
    L/2*norm(as.matrix(wL$beta-w$beta),type="F")^2+
    L/2*norm(as.matrix(wL$gamma-w$gamma),type="F")^2

  #while(criterion(w,treatment,mediators,outcome,mu=lambda/L,groups=groups,tau=tau)<criterion(wL,treatment,mediators,outcome,mu=lambda/L,groups=groups,tau=tau) & L<10000){
  while(test<0 & L<10000){

    L <-eta*L


    wgrad <- fgrad(w,treatment=treatment,mediators=mediators,outcome=outcome,tau=tau)
    u$alpha <- w$alpha - 1/L*wgrad$alpha
    u$beta <- w$beta - 1/L*wgrad$beta
    u$gamma <- w$gamma - 1/L*wgrad$gamma

    wL <- proximalapprox(u,mu=lambda/L,groups)

    test <- f(w,treatment,mediators,outcome,tau=tau)- f(wL,treatment,mediators,outcome,tau=tau) +
      sum(wgrad$alpha*(wL$alpha-w$alpha))+
      sum(wgrad$beta*(wL$beta-w$beta))+
      sum(wgrad$gamma*(wL$gamma-w$gamma))+
      L/2*norm((wL$alpha-w$alpha),type="F")^2+
      L/2*norm(as.matrix(wL$beta-w$beta),type="F")^2+
      L/2*norm(as.matrix(wL$gamma-w$gamma),type="F")^2

  }

  return(wL)

}
