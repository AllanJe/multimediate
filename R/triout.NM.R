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

