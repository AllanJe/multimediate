make.similarity = function(medata, similarity) {
  N <- ncol(medata)
  S <- matrix(rep(NA,N^2), ncol=N)
  for(i in 1:N) {
    for(j in 1:N) {
      S[i,j] <- similarity(medata[,i],medata[,j])
    }
  }
  return(S)
}
