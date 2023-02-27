#'@import stats
#'@importFrom graphics par
#'@importFrom graphics barplot
#'
SpecClust = function(medata,seuil=0){
  S <- make.similarity(medata, similarity)
  W <- make.affinity(S, seuil)
  D <- diag(apply(W, 1, sum))
  UL <- D - W
  round(UL,1)
  evL <- eigen(UL, symmetric=TRUE)
  par(mfrow=c(1,1))
  barplot(evL$values)
  k=as.numeric(readline(prompt="Nombre de cluster? "))
  Z   <- evL$vectors[,(ncol(evL$vectors)-k+1):ncol(evL$vectors)]
  return(kmeans(Z,k,nstart=3)$cluster)
}
