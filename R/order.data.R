order.data = function(Exposure,Outcome,Mediators,P){
  J=dim(Mediators)[1]#length(Exposure)
  K=dim(Mediators)[2]

  z=matrix(0,(K+1)*J,1)
  X=matrix(0,(K+1)*J,(P+1)*K+P)
  i=1
  p=1
  for(k in 1:K){
    #print(paste("Pour la colonne",k,"on remplit de",i,"à",(i+J-1),sep=' '))
    z[i:(i+J-1),1]=Mediators[,k]
    X[i:(i+J-1),p:(p+P-1)]=as.matrix(Exposure[1:J,1:P])
    i=i+J
    p=p+P
  }

  z[i:(i+J-1),1]=Outcome
  X[seq((K+1)*J-J+1,(K+1)*J),seq(P*K+1,(P+1)*K)]=as.matrix(Mediators)
  X[seq((K+1)*J-J+1,(K+1)*J),((P+1)*K+1):((P+1)*K+P)]=as.matrix(Exposure)

  group=c(rep(seq(1,K),each=P),seq(1,K),K+1)
  X=X[,order(group)] # on reordonne X pour gglasso
  group=c(rep(seq(1,K),each=P+1),K+1)
  return(list(z,X,group))
}

# order.data=function(Exposure,Outcome,Mediators){
#   J=length(Exposure)
#   K=dim(Mediators)[2]
#
#   z=matrix(0,(K+1)*J,1)
#   X=matrix(0,(K+1)*J,2*K+1)
#   i=1
#   for(k in 1:K){
#     #print(paste("Pour la colonne",k,"on remplit de",i,"à",(i+J-1),sep=' '))
#     z[i:(i+J-1),1]=Mediators[,k]
#     X[i:(i+J-1),k]=Exposure
#     i=i+J
#   }
#
#   z[i:(i+J-1),1]=Outcome
#   X[seq((K+1)*J-J+1,(K+1)*J),seq((K+1),K+K)]=as.matrix(Mediators)
#   X[seq((K+1)*J-J+1,(K+1)*J),2*K+1]=Exposure
#
#   group=c(rep(seq(1,K),2),K+1)
#   X=X[,order(group)] # on reordonne X pour gglasso
#   group=c(rep(seq(1,K),each=2),K+1)
#   return(list(z,X,group))
# }
