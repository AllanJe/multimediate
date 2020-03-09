triout.ci.NM <- function(d1.ci,p1.ci,d0.ci,p0.ci){
  NM=length(d1.ci[,1])
  res=array(NA,dim=c(4*NM,2))
  u=seq(1,4*NM,4)
  v=u+1
  w=u+2
  m=u+3

  for (i in 1:NM){
    res[u[i],] = d1.ci[i,]
    res[v[i],] = p1.ci[i,]
    res[w[i],] = d0.ci[i,]
    res[m[i],] = p0.ci[i,]
  }

  return(res)
}
