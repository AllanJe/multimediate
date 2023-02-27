triout.avg.NM = function(d,pm){
  NM=length(d)
  out=rep(NA,2*NM)
  u=seq(1,2*NM,2)
  out[u]=d
  out[u+1]=pm
  return(out)
}
