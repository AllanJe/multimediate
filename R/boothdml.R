boothdml <- function(treatment,mediators,outcome,lambda,Nboot=100,L0=.1,eta=2,tau=1,epsilon=.001){

  start=Sys.time()

  n <- dim(treatment)[1]
  P <- dim(treatment)[2]
  K <- dim(mediators)[2]

  if (n!=dim(mediators)[1]){
    print("The number of observations between mediators and treatment does not correspond")
  }

  if (n!=length(outcome)){
    print("The number of observations between outcome and treatment does not correspond")
  }



  counts <- rep(0,K)

  for(j in 1:Nboot){
    #construction des données bootstrap
    #print(j)
    bootsample <- sample(c(1:n),floor(n/2),replace=FALSE)
    boottreatment <- as.matrix(treatment[bootsample,])
    bootmediators <- as.matrix(mediators[bootsample,])
    bootoutcome <- outcome[bootsample]
    counts <- counts+ singlehdml(boottreatment,bootmediators,bootoutcome,groups=NULL,lambda=lambda,L0=L0,eta=eta,tau=tau,epsilon=epsilon)$selected
  }

  stop=Sys.time()
  return(list(counts=counts,time=stop-start))
}
