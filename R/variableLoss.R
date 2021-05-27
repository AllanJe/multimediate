variableLoss<- function(estimation,truth,binary=FALSE){ #returns the quadratic loss for non binary variables, the logistic one for binary variables
  
  return((1-binary)*sum((estimation-truth)^2)+binary*sum(-estimation*truth+log(1+exp(estimation))))
} 