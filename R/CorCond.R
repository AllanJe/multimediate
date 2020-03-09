#'
#'
#'@import stats
#'
#'
#'

CorCond=function(e,lmodel.m){
  NM=length(lmodel.m)
  k=1

  out=list()
  while(k < NM){
    for(l in (k+1):NM){
      a=paste("cor",names(lmodel.m[[k]]$model)[1],names(lmodel.m[[l]]$model)[1],sep=".")
      #print(a)
      out[[a]]=CorCond2(e,lmodel=list(lmodel.m[[k]],lmodel.m[[l]]))
    }
    k=k+1
  }
  Variance=rep(1,NM)
  for (z in 1:NM){
    if (inherits(lmodel.m[[z]],"lm") & !inherits(lmodel.m[[z]],"glm")){
      Variance[z]=var(lmodel.m[[z]]$residuals)
    }
  }
  sigmaestim=diag(Variance)

  i=1
  l=1
  vcor=c()
  while (i <NM){
    j=i+1
    for(k in j:NM){
      sigmaestim[i,k]=sigmaestim[k,i]=weighted.mean(out[[l]]$Correlation.estim,out[[l]]$Proportion,na.rm = TRUE)
      vcor=c(vcor,weighted.mean(out[[l]]$Correlation.estim,out[[l]]$Proportion,na.rm = TRUE))
      l=l+1
    }
    i=i+1
  }
  out[["sigmaestim"]]=sigmaestim
  out[["vcor"]]=vcor
  return(out)
}
