pval = function(x, xhat,seu=0){
  if (xhat == seu) {out <- 1}
  else {
    out <- 2 * min(sum(x > seu,na.rm=TRUE), sum(x < seu,na.rm=TRUE)) / length(x)
  }
  return(min(out, 1))
}
