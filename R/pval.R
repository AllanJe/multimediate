pval <- function(x, xhat){
  if (xhat == 0) {out <- 1}
  else {
    out <- 2 * min(sum(x > 0,na.rm=TRUE), sum(x < 0,na.rm=TRUE)) / length(x)
  }
  return(min(out, 1))
}
