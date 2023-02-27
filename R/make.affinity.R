make.affinity = function(S, seuil) {
  N <- length(S[,1])

  if (seuil==0) {  # fully connected
    W <- S
  }

  else
  {
    W <- matrix(rep(0,N^2), ncol=N)
    for(i in 1:N)
    { # for each line
      # only connect to those points with larger similarity
      #best.similarities <- sort(S[i,], decreasing=TRUE)
      #for (s in best.similarities) {
      for(j in 1:N) { #j <- which(S[i,] == s)
        if( S[i,j]>=seuil){
          W[i,j] <- S[i,j]
          W[j,i] <- S[i,j] # to make an undirected graph, ie, the matrix becomes symmetric
        }
      }
    }
  }
  return(W)
}
