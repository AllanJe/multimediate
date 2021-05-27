#' lambdaChoice
#'
#' "lambdaChoice" is used to .
#'
#'
#' @param treatment  qsdfgh
#' @param mediators ertyu
#' @param outcome  poiu
#' @param lambdamax mlkjh
#' @param N1 kjhg
#' @param selectedMin mlkjh
#' @param selectedMax zerf
#' @param L0 zert
#' @param eta dh
#' @param tau sdfghj
#' @param epsilon cdfg
#'
#' @return lambda value
#' @export
#'

lambdaChoice <- function(treatment,mediators,outcome,lambdamax=500,N1=10,selectedMin=.5*length(outcome),selectedMax=length(outcome),L0=1,eta=2,tau=1,epsilon=.001){

  lambdamin=0
  lambda=lambdamax


  selected <- c()
  for (i in 1:N1){
    selected <- c(selected,singlehdml(treatment,mediators,outcome,groups=NULL,lambda=lambda,L0=L0,eta=eta,tau=tau,epsilon=epsilon))$selectedNumber
  }



  if(mean(selected)>selectedMax){
    print("choose a greater lambdamax")
  }

  else{

    i=1       # ne pas faire plus plus de 10 divisions
    lambdamax=lambda
    test=FALSE
    while((test==FALSE) & (i<10)){

      i=i+1
      lambda=mean(c(lambdamin,lambdamax))
      selected <- c()
      for (j in 1:N1){
        selected <- c(selected,singlehdml(treatment,mediators,outcome,groups=NULL,lambda=lambda,L0=L0,eta=eta,tau=tau,epsilon=epsilon))$selectedNumber
      }

      if(mean(selected)<selectedMin){ lambdamax <- lambda}
      else if(mean(selected)>selectedMax){ lambdamin <- lambda}
      else{test <- TRUE}
    }


  }

  return(lambda)

}
