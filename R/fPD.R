#' @importFrom mvtnorm pmvnorm
fPD = function(rho,covMiMj,pMi,bornep,p,pMj,borneq,q){
  Sigma1=matrix(c(1,rho,rho,1),2,2)
  pmvn=NULL
  for (i in 1:length(p)){
    for(j in 1:length(q)){
      pmvn=c(pmvn,p[i]*q[j]*pmvnorm(lower=c(bornep[i],borneq[j]),upper=c(bornep[i+1],borneq[j+1]),sigma=Sigma1))
    }
  }
  out=covMiMj-sum(pmvn)+pMi*pMj
  return(out)
}
