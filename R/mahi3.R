#' mahi3
#'
#'
#' 'mahi3' is used to select mediators trough a high number of mediators and estimate various quantities for causal mediation analysis using "multimediate", with the selected mediators.
#'
#'@param data a data.frame with the exposure, the outcome and mediators.
#'@param name.exposure a character string indicating the name of the exposure in data.
#'@param name.outcome a character string indicating the name of the outcome in data.
#'@param name.mediators a character string indicating the name of mediators in data.
#'@param Nboot number of bootstrap for stability selection of mediators. Default is 100.
#'@param lambda vector or a single value of lambda use for group-lasso.
#'@param L0 dfghj
#'@param eta dfghj
#'@param tau dfghj
#'@param epsilon cvbn
#'@param Kmax maximum of mediators keept after the ranking obtained with stability selection.
#'@param p.adjust.method a character string indicating the name of the method for the correction for the multiple test. See help with p.adjust.methods.
#'@param pvalseuil the p-value for the multiple test
#'@param bin a logical value indicating if the outcome is binary. if 'TRUE' a probit regression will be use is the second step. Default is 'FALSE'.
#'@param step2 a logical value indicating if the a multiple test using the p-value calculated with multimediate (a multiple mediation analysis) have to be made. Default is 'FALSE'.
#'@return mahi returns an object of class "mahi", a list that contains the components listed below.
#'
#'
#' @export
#'

mahi3=function(data,name.exposure,name.outcome,name.mediators,lambda,Nboot=100,L0=.1,eta=2,tau=1/P,epsilon=.001,Kmax=NULL,bin=FALSE,step2=FALSE,p.adjust.method="bonferroni",pvalseuil=0.05){
  n=dim(data)[1]
  K=length(name.mediators)
  P=length(name.exposure)
  if(is.null(Kmax)){
    Kmax=2*n/log(n)
  }


  treatment=as.matrix(data[,name.exposure])
  mediators=as.matrix(data[,name.mediators])
  outcome=data[,name.outcome]

  bootstep <- boothdml(treatment,mediators,outcome,lambda,Nboot,L0,eta,tau,epsilon)
  ranking.id=order(bootstep$counts,decreasing = TRUE)
  ranking=name.mediators[ranking.id]
  Kmaxrank= ranking[1:Kmax]

if(!step2){
  out=list(Kmaxrank=Kmaxrank,ranking=ranking,ranking.id=ranking.id,bootcount=bootstep$counts,
           n=n,K=K,Kmax=Kmax,lambda=lambda)
}
else{ print("Step 2: Multiple test on indirect effects")
  selection=rep(FALSE,K)
  if(P>1){
    selectionbyt=array(FALSE,dim=c(K,P))

  }
  multimed =pvals=pvalscorr=list()
  pvalsbytreat=matrix(NA,Kmax,P)
  for(p in 1:P){
    lmodel.m=list()
    medforout=""
    for(medch in Kmaxrank){
      lmodel.m[[medch]] = lm(formula(paste( medch,"~",name.exposure[p])), data = data)
      medforout=paste(medforout,medch,sep=" + ")
    }
    formout1 = paste(name.outcome,"~",name.exposure[p])
    if(bin){
      model.y = glm(formula( paste(formout1,medforout,sep=" ")), data = data, family=binomial(link = "logit"))
    }
    else{model.y = lm(formula( paste(formout1,medforout,sep=" ")), data = data)}


    correlated=TRUE
    if(length(Kmaxrank)==1){correlated=FALSE}
    multimed[[p]]=multimediate(lmodel.m,correlated,model.y,treat=name.exposure[p])
    pvalsbytreat[,p]=summary(multimed[[p]],opt="avg")[seq(3,length(Kmaxrank)*2+2,2),5]
  }

  pvalscorr=apply(pvalsbytreat,2,p.adjust,method = p.adjust.method)
  selmed1=pvalscorr<=pvalseuil
  if(P>1){
    for (p in 1:P){
      selmedbyt=Kmaxrank[selmed1[,p]]
      selectionbyt[which(name.mediators %in% selmedbyt),p]=TRUE
    }}

  selmed2=apply(selmed1,1,sum)
  selmed3=Kmaxrank[selmed2==P]
  selection[which(name.mediators %in% selmed3)]=TRUE




  if(P>1){
    out=list(selectionbyt=selectionbyt,selection=selection,
             Kmaxrank=Kmaxrank,ranking=ranking,ranking.id=ranking.id,bootcount=bootstep$counts,
             multimed=multimed,pvals=pvalsbytreat,pvalscorr=pvalscorr,
             n=n,K=K,Kmax=Kmax,lambda=lambda)
  }
  else{
    out=list(selection=selection,
             Kmaxrank=Kmaxrank,ranking=ranking,ranking.id=ranking.id,bootcount=bootstep$counts,
             multimed=multimed[[1]],pvals=pvalsbytreat,pvalscorr=pvalscorr,
             n=n,K=K,Kmax=Kmax,lambda=lambda)
  }}
}
