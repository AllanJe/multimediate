triout.ci.avg.NM <- function(d.ci,pm.ci){
  NM=length(d.ci[,1])
  res=array(NA,dim=c(2*NM,2))
  u=seq(1,2*NM,2)
  v=u+1

  for (i in 1:NM){
    res[u[i],] = d.ci[i,]
    res[v[i],] = pm.ci[i,]
  }

  return(res)
}
