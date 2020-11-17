#' mahi
#'
#'
#' 'mahi' is used to select mediators trough a high number of mediators and estimate various quantities for causal mediation analysis using "multimediate", with the selected mediators.
#'
#'@param data a data.frame with the exposure, the outcome and mediators.
#'@param name.exposure a character string indicating the name of the exposure in data.
#'@param name.outcome a character string indicating the name of the outcome in data.
#'@param name.mediators a character string indicating the name of mediators in data.
#'@param nboot number of bootstrap for stability selection of mediators. Default is 100.
#'@param nlambda number of lambda use when lambda is not specified.
#'@param lambda vector or a single value of lambda use for group-lasso.
#'@param Kmax maximum of mediators keept after the ranking obtained with stability selection.
#'@param p.adjust.method a character string indicating the name of the method for the correction for the multiple test. See help with p.adjust.methods.
#'@param pvalseuil the p-value for the multiple test
#@param bin a logical value indicating if the outcome is binary. if 'TRUE' a probit regression will be use is the second step. Default is 'FALSE'.
#'@return mahi returns an object of class "mahi", a list that contains the components listed below.
#'
#' @export
#'
#' @importFrom tsutils lambdaseq
#' @importFrom gglasso gglasso

mahi=function(data,name.exposure,name.outcome,name.mediators,nboot=100,nlambda=50,lambda=NULL,Kmax=NULL,p.adjust.method="bonferroni",pvalseuil=0.05){#,bin=FALSE){
  n=dim(data)[1]
  K=length(name.mediators)
  P=length(name.exposure)

  if(!is.null(lambda)){
    nlambda <- length(lambda)
  }

  selectMed=array(NA,dim=c(K,nlambda,nboot))

  if(is.null(Kmax)){
    Kmax=2*n/log(n)
  }


  if(is.null(lambda)){

    if (P==1){
      Exposure=matrix(data[,name.exposure],n,1)
    }
    else {
      Exposure=data[,name.exposure]
    }
    Outcome=data[,name.outcome]
    Mediators=data[,name.mediators]
    test=order.data(Exposure,Outcome,Mediators,P)
    z=test[[1]]
    X=test[[2]]
    lambda <- lambdaseq(X,z,nLambda=nlambda)$lambda
  }
  print("Etape 1: GG-LASSO")
  for(bt in 1:nboot){
    T1=Sys.time()
    print(paste(bt,"/",nboot))
    databt=data[sample(1:n,n,replace = TRUE),]
    if (P==1){
      Exposure=matrix(databt[,name.exposure],n,1)
    }
    else {
      Exposure=databt[,name.exposure]
    }
    Outcome=databt[,name.outcome]
    Mediators=databt[,name.mediators]
    test=order.data(Exposure,Outcome,Mediators,P)
    z=test[[1]]
    X=test[[2]]
    group=test[[3]]
    print("GG-LASSO")
    ggltest=gglasso(x=X,y=z,group=group,lambda=lambda)
    print("Select")
    for(k in 1:K){
      u=seq(2,(P+1)*K,P+1)[k]
      selectMed[k,,bt]=(coef(ggltest)[u,]!=0)# & coef(ggltest)[u+1,]!=0)
    }
    T2=Sys.time()
    Times=difftime(T2,T1)
    print(paste("L'etape 1 sera terminee vers :",T1+Times*(nboot-bt+1),sep=""))
  }


  rang=apply(selectMed,c(1,2),sum)
  etap1=matrix(NA,Kmax,nlambda)
  for(lamb in 1:nlambda){
    etap1[,lamb]=order(rang[,lamb],decreasing = TRUE)[1:Kmax]
  }


  print("Etape 2: Tests multiples")
  selection=array(FALSE,dim=c(K,nlambda))
  if(P>1){
    selectionbyt=array(FALSE,dim=c(K,P,nlambda))
  }
  multimed =pvals=pvalscorr=pvalsbytreat=list()

  for(lamb in 1:nlambda){
    T1=Sys.time()
    print(paste(lamb,"/",nlambda))
    mediators = name.mediators[etap1[,lamb]]
    pvalsbytreat[[lamb]]=matrix(NA,length(mediators),P)
    multimed[[lamb]]=list()
    for(p in 1:P){
      lmodel.m=list()
      medforout=""
      for(medch in mediators){
        lmodel.m[[medch]] = lm(formula(paste( medch,"~",name.exposure[p])), data = data)
        medforout=paste(medforout,medch,sep=" + ")
      }
      formout1 = paste(name.outcome,"~",name.exposure[p])
      model.y = lm(formula( paste(formout1,medforout,sep=" ")), data = data)
      # if(bin==FALSE){
      #   model.y = lm(formula( paste(formout1,medforout,sep=" ")), data = data)
      # }
      # else {
      #   model.y = glm(formula( paste(formout1,medforout,sep=" ")), data = data,family = binomial("probit"))
      # }

      correlated=TRUE
      if(length(mediators)==1){correlated=FALSE}
      multimed[[lamb]][[p]]=multimediate(lmodel.m,correlated,model.y,treat=name.exposure[p])
      pvalsbytreat[[lamb]][,p]=summary(multimed[[lamb]][[p]],opt="avg")[seq(3,length(mediators)*2+2,2),5]
    }
    #pvals[[lamb]]=pvalsbytreat
    pvalscorr[[lamb]]=apply(pvalsbytreat[[lamb]],2,p.adjust,method = p.adjust.method)#p.adjust(pvals[[lamb]],method = p.adjust.method)
    selmed1=pvalscorr[[lamb]]<=pvalseuil
    if(P>1){
      for (p in 1:P){
        selmedbyt=mediators[selmed1[,p]]
        selectionbyt[which(name.mediators %in% selmedbyt),p,lamb]=TRUE
      }}

    selmed2=apply(selmed1,1,sum)
    selmed3=mediators[selmed2==P]
    selection[which(name.mediators %in% selmed3),lamb]=TRUE

    T2=Sys.time()
    Times=difftime(T2,T1)
    print(paste("L'etape 2 sera terminee avant :",T1+Times*(nlambda-lamb+1),sep=""))
  }


  if(P>1){
    out=list(multimed=multimed,selectMed=selectMed,selection=selection,selectionbyt=selectionbyt,n=n,K=K,Kmax=Kmax,lambda=lambda,nlambda=nlambda,
             rang=rang,etap1=etap1,pvals=pvals,pvalscorr=pvalscorr)
  }
  else{
    out=list(multimed=multimed,selectMed=selectMed,selection=selection,n=n,K=K,Kmax=Kmax,lambda=lambda,nlambda=nlambda,
             rang=rang,etap1=etap1,pvals=pvals,pvalscorr=pvalscorr)
  }

}

# mahi=function(data,name.exposure,name.outcome,name.mediators,nboot=30,nlambda=1,lambda=NULL,Kmax=NULL,p.adjust.method="bonferroni",pvalseuil=0.05){
#   n=dim(data)[1]
#   K=length(name.mediators)
#
#   if(!is.null(lambda)){
#     nlambda <- length(lambda)
#   }
#
#   selectMed=array(NA,dim=c(K,nlambda,nboot))
#
#   if(is.null(Kmax)){
#     Kmax=2*n/log(n)
#   }
#
#
#   if(is.null(lambda)){
#     Exposure=data[,name.exposure]
#     Outcome=data[,name.outcome]
#     Mediators=data[,name.mediators]
#     test=order.data(Exposure,Outcome,Mediators)
#     z=test[[1]]
#     X=test[[2]]
#     lambda <- lambdaseq(X,z,nLambda=nlambda)$lambda
#   }
#   print("Etape 1: GG-LASSO")
#   for(bt in 1:nboot){
#     T1=Sys.time()
#     print(paste(bt,"/",nboot))
#     databt=data[sample(1:n,n,replace = TRUE),]
#     Exposure=databt[,name.exposure]
#     Outcome=databt[,name.outcome]
#     Mediators=databt[,name.mediators]
#     test=order.data(Exposure,Outcome,Mediators)
#     z=test[[1]]
#     X=test[[2]]
#     group=test[[3]]
#
#     ggltest=gglasso(x=X,y=z,group=group,lambda=lambda)
#
#     for(k in 1:K){
#       u=seq(2,2*K,2)[k]
#       selectMed[k,,bt]=(coef(ggltest)[u,]!=0 & coef(ggltest)[u+1,]!=0)
#     }
#     T2=Sys.time()
#     Times=difftime(T2,T1)
#     print(paste("L'etape 1 sera terminee vers :",T1+Times*(nboot-bt+1),sep=""))
#   }
#
#
#   rang=apply(selectMed,c(1,2),sum)
#   etap1=matrix(NA,Kmax,nlambda)
#   for(lamb in 1:nlambda){
#     etap1[,lamb]=order(rang[,lamb],decreasing = TRUE)[1:Kmax]
#   }
#
#
#   print("Etape 2: Tests multiples")
#   selection=array(FALSE,dim=c(K,nlambda))
#   multimed =pvals=pvalscorr=list()
#   for(lamb in 1:nlambda){
#     T1=Sys.time()
#     print(paste(lamb,"/",nlambda))
#     mediators = name.mediators[etap1[,lamb]]
#     lmodel.m=list()
#     medforout=""
#     for(medch in mediators){
#       lmodel.m[[medch]] = lm(formula(paste( medch,"~",name.exposure)), data = data)
#       medforout=paste(medforout,medch,sep=" + ")
#     }
#     formout1 = paste(name.outcome,"~",name.exposure)
#     model.y = lm(formula( paste(formout1,medforout,sep=" ")), data = data)
#
#     correlated=TRUE
#     if(length(mediators)==1){correlated=FALSE}
#     multimed[[lamb]]=multimediate(lmodel.m,correlated,model.y,treat=name.exposure)
#
#     pvals[[lamb]]=summary(multimed[[lamb]],opt="avg")[seq(3,length(mediators)*2+2,2),5]
#     pvalscorr[[lamb]]=p.adjust(pvals[[lamb]],method = p.adjust.method)
#     selmed=mediators[pvalscorr[[lamb]]<=pvalseuil]
#     selection[which(name.mediators %in% selmed),lamb]=TRUE
#
#     T2=Sys.time()
#     Times=difftime(T2,T1)
#     print(paste("L'etape 2 sera terminee avant :",T1+Times*(nlambda-lamb+1),sep=""))
#   }
#
#   out=list(multimed=multimed,selectMed=selectMed,selection=selection,n=n,K=K,Kmax=Kmax,lambda=lambda,nlambda=nlambda,
#            rang=rang,etap1=etap1,pvals=pvals,pvalscorr=pvalscorr)
#
# }
