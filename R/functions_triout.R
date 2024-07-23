triout.NM = function(d1,pm1,d0,pm0){
  NM=length(d1)
  out=rep(NA,4*NM)
  u=seq(1,4*NM,4)
  out[u]=d1
  out[u+1]=pm1
  out[u+2]=d0
  out[u+3]=pm0
  return(out)
}

triout.ci.NM = function(d1.ci,p1.ci,d0.ci,p0.ci){
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

triout.ci.avg.NM = function(d.ci,pm.ci){
  NM = length(d.ci[,1])
  res = array(NA,dim=c(2*NM,2))
  u = seq(1,2*NM,2)
  v = u+1

  for (i in 1:NM){
    res[u[i],] = d.ci[i,]
    res[v[i],] = pm.ci[i,]
  }

  return(res)
}

triout.avg.NM = function(d,pm){
  NM=length(d)
  out=rep(NA,2*NM)
  u=seq(1,2*NM,2)
  out[u]=d
  out[u+1]=pm
  return(out)
}
