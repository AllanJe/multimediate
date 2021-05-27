singlehdml <- function(treatment,mediators,outcome,groups=c(),lambda,L0=1,eta=2,tau=1,epsilon=.001){

  n <- dim(treatment)[1]
  P <- dim(treatment)[2]
  K <- dim(mediators)[2]

  if (n!=dim(mediators)[1]){
    print("The number of observations between mediators and treatment does not correspond")
  }

  if (n!=length(outcome)){
    print("The number of observations between outcome and treatment does not correspond")
  }

  # ajout d'une colonne de 1 au début de T qui permet de prendre en compte les ordonnées à l'origine dans les deux combinaisons linéaires
  treatment <- cbind(rep(1,n),treatment)

  # initialisation de alpha, beta, gamma
  w <- c()
  w$alpha <- matrix(runif(K*P+K,-1,1),nrow = P+1,ncol=K)   #la première ligne contient les ordonnées à l'origine pour les Mk, \alpha_{pk} correspond à alpha[p+1,k]
  w$beta <- runif(K,-1,1)
  w$gamma <- runif(P+1,-1,1)               #la première ligne contient l'ordonnée à l'origine pour Y, \gamma_{p} correspond à gamma[p+1]



  # vecteurs des groupes pour retrouver les coefficients à prendre en compte dans la pénalité
  # pour alpha, la k^e colonne est dans le groupe k sauf l'ordonnée à l'origine qui est dans 0, pour beta, k^e coord dans k, tout gamma dans 0
  groups <- c()
  groups$alpha <- rbind(rep(0,K),matrix(c(1:K),nrow=P,ncol=K,byrow=TRUE))
  groups$beta <- c(1:K)
  groups$gamma <- rep(0,P+1)

  #print("start")
  #print(criterion(w,treatment,mediators,outcome,mu=lambda,groups=groups,tau=tau))
  wnew <- oneIteration(w,treatment,mediators,outcome,groups=groups,lambda=lambda,L0=L0,eta=eta,tau=tau)

  i=1
  while(criterion(w,treatment,mediators,outcome,mu=lambda,groups=groups,tau=tau)-criterion(wnew,treatment,mediators,outcome,mu=lambda,groups=groups,tau=tau)>epsilon & i<10){
    i=i+1
    #print(i)
    #print(criterion(w,treatment,mediators,outcome,mu=lambda,groups=groups,tau=tau))
    w <- wnew
    wnew <- oneIteration(w,treatment,mediators,outcome,groups=groups,lambda=lambda,L0=L0,eta=eta,tau=tau)
  }

  selectedNumber <- 0
  selected <- rep(0,K)
  for(k in 1:K){
    knorm=norm(as.matrix( c(wnew$alpha[groups$alpha==k],wnew$beta[groups$beta==k])),type="F")
    selected[k] <- (knorm>0)*1
    selectedNumber <- selectedNumber + (knorm>0)
  }

  return(list(wnew=wnew,selectedNumber=selectedNumber,selected=selected,criterion=criterion(wnew,treatment,mediators,outcome,mu=lambda,groups,tau=tau)))
}
