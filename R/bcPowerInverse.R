bcPowerInverse=function(U,lambda,gamma=NULL){

  U[U<(-1/lambda)]=NaN

  out=(lambda * U + 1)^(1/lambda)

  if(!is.null(gamma)){
    out= ( (4 * out^2) - gamma^2)/(4 * out)
  }

    return(out)
  }
