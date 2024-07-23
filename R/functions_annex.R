#' getvarnames
#'
#'
#' 'getvarnames' is used to get all variable names from a regression model
#'
#'@param formula the formula from which to extract variable names
#'@param data data frame with variables from the formula
#'@return 'getvarnames' returns a list with 'varnames' (referring to all variable names), 'xvar' (referring to the predictors), and 'yvar' (referring to the outcome)

getvarnames = function(formula, data = NULL)
{
  if (is.character(formula))
    return(list(varnames=formula, xvar=formula, yvar=NULL))
  if (is.null(formula)) return(list(varnames=NULL, xvar=NULL, yvar=NULL))

  formula <- formula(formula)
  lyv <- NULL
  lxv <- lvnm <- all.vars(formula[1:2])
  if (length(formula)==3) {
    lyv <- lxv
    lxv <- all.vars(formula[-2])
    if ("." %in% lxv) {
      if (length(data)==0)
        stop("!getvarnames! '.' in formula and no 'data'")
      lform <- formula(terms(formula, data=data))
      lxv <- all.vars(lform[-2])
    }
    lvnm <- c(lxv, lvnm)
  }
  return(list(varnames=lvnm, xvar=lxv, yvar=lyv))
}

plot.process = function(model,logit="effects"){

  if (logit=="OR"){
    coef.vec.1 <- c(model$ORd1, model$ORz1)
    lower.vec.1 <- c(model$ORd1.ci[1], model$ORz1.ci[1])
    upper.vec.1 <- c(model$ORd1.ci[2], model$ORz1.ci[2])
    coef.vec.1.NM <- model$ORd1.NM
    lower.vec.1.NM <- model$ORd1.ci.NM[,1]
    upper.vec.1.NM <- model$ORd1.ci.NM[,2]
    tau.vec <- c(model$ORtau.coef, model$ORtau.ci[1], model$ORtau.ci[2])
    range.1 <- range(model$ORd1.ci[1],model$ORd1.ci.NM[,1], model$ORz1.ci[1], model$ORtau.ci[1],
                     model$ORd1.ci[2],model$ORd1.ci.NM[,2], model$ORz1.ci[2], model$ORtau.ci[2])

    coef.vec.0 <- c(model$ORd0, model$ORz0)
    lower.vec.0 <- c(model$ORd0.ci[1], model$ORz0.ci[1])
    upper.vec.0 <- c(model$ORd0.ci[2], model$ORz0.ci[2])
    coef.vec.0.NM <- model$ORd0.NM
    lower.vec.0.NM <- model$ORd0.ci.NM[,1]
    upper.vec.0.NM <- model$ORd0.ci.NM[,2]
    range.0 <- range(model$ORd0.ci[1],model$ORd0.ci.NM[,1], model$ORz0.ci[1], model$ORtau.ci[1],
                     model$ORd0.ci[2],model$ORd0.ci.NM[,2], model$ORz0.ci[2], model$ORtau.ci[2])

    coef.vec.avg <- c(model$ORd.avg, model$ORz.avg)
    lower.vec.avg <- c(model$ORd.avg.ci[1], model$ORz.avg.ci[1])
    upper.vec.avg <- c(model$ORd.avg.ci[2], model$ORz.avg.ci[2])
    coef.vec.avg.NM <- model$ORd.avg.NM
    lower.vec.avg.NM <- model$ORd.avg.ci.NM[,1]
    upper.vec.avg.NM <- model$ORd.avg.ci.NM[,2]
    range.avg <- range(model$ORd.avg.ci[1], model$ORd.avg.ci.NM[,1],model$ORz.avg.ci[1], model$ORtau.ci[1],
                       model$ORd.avg.ci[2],model$ORd.avg.ci.NM[,2], model$ORz.avg.ci[2], model$ORtau.ci[2])
  }
  else if (logit=="logOR"){
    coef.vec.1 <- c(model$logORd1, model$logORz1)
    lower.vec.1 <- c(model$logORd1.ci[1], model$logORz1.ci[1])
    upper.vec.1 <- c(model$logORd1.ci[2], model$logORz1.ci[2])
    coef.vec.1.NM <- model$logORd1.NM
    lower.vec.1.NM <- model$logORd1.ci.NM[,1]
    upper.vec.1.NM <- model$logORd1.ci.NM[,2]
    tau.vec <- c(model$logORtau.coef, model$logORtau.ci[1], model$logORtau.ci[2])
    range.1 <- range(model$logORd1.ci[1],model$logORd1.ci.NM[,1], model$logORz1.ci[1], model$logORtau.ci[1],
                     model$logORd1.ci[2],model$logORd1.ci.NM[,2], model$logORz1.ci[2], model$logORtau.ci[2])

    coef.vec.0 <- c(model$logORd0, model$logORz0)
    lower.vec.0 <- c(model$logORd0.ci[1], model$logORz0.ci[1])
    upper.vec.0 <- c(model$logORd0.ci[2], model$logORz0.ci[2])
    coef.vec.0.NM <- model$logORd0.NM
    lower.vec.0.NM <- model$logORd0.ci.NM[,1]
    upper.vec.0.NM <- model$logORd0.ci.NM[,2]
    range.0 <- range(model$logORd0.ci[1],model$logORd0.ci.NM[,1], model$logORz0.ci[1], model$logORtau.ci[1],
                     model$logORd0.ci[2],model$logORd0.ci.NM[,2], model$logORz0.ci[2], model$logORtau.ci[2])

    coef.vec.avg <- c(model$logORd.avg, model$logORz.avg)
    lower.vec.avg <- c(model$logORd.avg.ci[1], model$logORz.avg.ci[1])
    upper.vec.avg <- c(model$logORd.avg.ci[2], model$logORz.avg.ci[2])
    coef.vec.avg.NM <- model$logORd.avg.NM
    lower.vec.avg.NM <- model$logORd.avg.ci.NM[,1]
    upper.vec.avg.NM <- model$logORd.avg.ci.NM[,2]
    range.avg <- range(model$logORd.avg.ci[1], model$logORd.avg.ci.NM[,1],model$logORz.avg.ci[1], model$logORtau.ci[1],
                       model$logORd.avg.ci[2],model$logORd.avg.ci.NM[,2], model$logORz.avg.ci[2], model$logORtau.ci[2])
  }
  else{
    coef.vec.1 <- c(model$d1, model$z1)
    lower.vec.1 <- c(model$d1.ci[1], model$z1.ci[1])
    upper.vec.1 <- c(model$d1.ci[2], model$z1.ci[2])
    coef.vec.1.NM <- model$d1.NM
    lower.vec.1.NM <- model$d1.ci.NM[,1]
    upper.vec.1.NM <- model$d1.ci.NM[,2]
    tau.vec <- c(model$tau.coef, model$tau.ci[1], model$tau.ci[2])
    range.1 <- range(model$d1.ci[1],model$d1.ci.NM[,1], model$z1.ci[1], model$tau.ci[1],
                     model$d1.ci[2],model$d1.ci.NM[,2], model$z1.ci[2], model$tau.ci[2])

    coef.vec.0 <- c(model$d0, model$z0)
    lower.vec.0 <- c(model$d0.ci[1], model$z0.ci[1])
    upper.vec.0 <- c(model$d0.ci[2], model$z0.ci[2])
    coef.vec.0.NM <- model$d0.NM
    lower.vec.0.NM <- model$d0.ci.NM[,1]
    upper.vec.0.NM <- model$d0.ci.NM[,2]
    range.0 <- range(model$d0.ci[1],model$d0.ci.NM[,1], model$z0.ci[1], model$tau.ci[1],
                     model$d0.ci[2],model$d0.ci.NM[,2], model$z0.ci[2], model$tau.ci[2])

    coef.vec.avg <- c(model$d.avg, model$z.avg)
    lower.vec.avg <- c(model$d.avg.ci[1], model$z.avg.ci[1])
    upper.vec.avg <- c(model$d.avg.ci[2], model$z.avg.ci[2])
    coef.vec.avg.NM <- model$d.avg.NM
    lower.vec.avg.NM <- model$d.avg.ci.NM[,1]
    upper.vec.avg.NM <- model$d.avg.ci.NM[,2]
    range.avg <- range(model$d.avg.ci[1], model$d.avg.ci.NM[,1],model$z.avg.ci[1], model$tau.ci[1],
                       model$d.avg.ci[2],model$d.avg.ci.NM[,2], model$z.avg.ci[2], model$tau.ci[2])
  }






  return(list(coef.vec.1 = coef.vec.1, lower.vec.1 = lower.vec.1, upper.vec.1 = upper.vec.1,coef.vec.1.NM=coef.vec.1.NM,lower.vec.1.NM=lower.vec.1.NM, upper.vec.1.NM= upper.vec.1.NM,
              coef.vec.0 = coef.vec.0, lower.vec.0 = lower.vec.0, upper.vec.0 = upper.vec.0,coef.vec.0.NM=coef.vec.0.NM,lower.vec.0.NM=lower.vec.0.NM, upper.vec.0.NM= upper.vec.0.NM,
              coef.vec.avg = coef.vec.avg, lower.vec.avg = lower.vec.avg, upper.vec.avg = upper.vec.avg,coef.vec.avg.NM=coef.vec.avg.NM,lower.vec.avg.NM=lower.vec.avg.NM, upper.vec.avg.NM= upper.vec.avg.NM,
              tau.vec = tau.vec,
              range.1 = range.1, range.0 = range.0, range.avg = range.avg))
}

pval = function(x, xhat,seu=0){
  if (is.na(xhat)) {
    out =NA
  }
  else{
    if (xhat == seu) {out <- 1}
    else {
      out <- 2 * min(sum(x > seu,na.rm=TRUE), sum(x < seu,na.rm=TRUE)) / length(x)
    }
  }

  return(min(out, 1))

}


#' @import stats
#' @importFrom MASS mvrnorm

sim_multi=function(n=1000, NM=3,
                   binaryOutcome = FALSE,
                   coeff, correlation=c(0.9,0.6,0.3),
                   covariables){
  Treatment=sample(0:1,n,replace=TRUE)
  data=data.frame(Intercept=1,Treatment=Treatment)


  lengthcovariables <- dim(covariables)[1]
  if(!is.null(covariables)){
    for(lc in 1:dim(covariables)[2]){
      if( length(covariables[,lc])!= n ){
        stop("Not all covariables are of length n.")
      }
      assign(names(covariables)[lc],covariables[,lc])
    }
    data=cbind(data,covariables)

  }



  nuplet=unique(data)



  ntre=c("Intercept","Treatment")
  nmed=paste("M",seq(1,NM),sep="")
  row.names(coeff)=c(nmed,"Outcome")

  if(!is.null(covariables)){
    ncov=paste("C",seq(1,dim(covariables)[2]),sep="")
    colnames(coeff)=c(ntre,ncov,nmed)
  }
  else{
    colnames(coeff)=c(ntre,nmed)
  }




  Valuev=NULL
  for(i in 1:dim(nuplet)[1]){
    Valuet=NULL
    for (j in 1:(dim(nuplet)[2]-1)){
      Valuet=paste(Valuet,paste(nuplet[i,-1])[j],sep="")
    }
    Valuev=c(Valuev,Valuet)
  }
  cf.med.names=NULL
  for(nm in 1:NM){
    cf.med.names=c(cf.med.names,paste(nmed[nm],Valuev,sep="."))
  }

  muM=NULL
  for (nm in 1:NM){
    for(i in 1:dim(nuplet)[1]){
      muM=c(muM,sum(coeff[nm,1:(ncol(coeff)-NM)]*nuplet[i,]))
    }
  }

  Sigma=matrix(0,length(muM),length(muM))
  u=nm=1
  v=dim(nuplet)[1]
  for(nm in 1:NM){
    Sigma[u:v,u:v]=matrix(1,dim(nuplet)[1],dim(nuplet)[1])
    nm=nm+1
    u=v+1
    v=dim(nuplet)[1]*nm
  }

  a=seq(1,length(muM),by=dim(nuplet)[1])

  l=1
  i=1
  while(i<NM){
    u=a[i]
    for(k in (i+1):(length(a))){
      v=a[k]
      w=v+dim(nuplet)[1]
      Sigma[u:(u+dim(nuplet)[1]-1),v:(w-1)]=Sigma[v:(w-1),u:(u+dim(nuplet)[1]-1)]=matrix(correlation[l],dim(nuplet)[1],dim(nuplet)[1])
      l=l+1
    }
    i=i+1
  }

  cf.data=MASS::mvrnorm(n=n,muM,Sigma,tol=10)
  cf.data=as.data.frame(cf.data)
  names(cf.data)=cf.med.names
  cf.data$Treatment=Treatment



  for (i in 1:dim(nuplet)[1]){
    transct=cf.data[,paste(nmed,Valuev[i],sep=".")]
    if (binaryOutcome==TRUE){
      Yerror <- stats::rlogis(n,0,1)
      cf.data[,paste("Y",Valuev[i],sep=".")]=apply(coeff[NM+1,]*t(cbind(nuplet[i,],transct)),2,sum)+Yerror>0
    }
    else{
      Yerror <- stats::rnorm(n,0,1)
      cf.data[,paste("Y",Valuev[i],sep=".")]=apply(coeff[NM+1,]*t(cbind(nuplet[i,],transct)),2,sum)+Yerror
    }

  }


  for(nm in 1:NM){
    data[nmed[nm]]=NA
  }
  data$Outcome=NA
  data=data[,-1]

  com=c()
  for (i in 1:dim(nuplet)[1]){
    com = c(com,parse(text= paste(paste(names(nuplet)[-1],nuplet[i,-1], sep = "=="), collapse = " & ")))
  }
  for (i in 1:dim(nuplet)[1]){
    for(nm in 1:NM){
      data[nmed[nm]][which(eval(com[i])),]=cf.data[,paste(nmed[nm],Valuev[i],sep=".")][which(eval(com[i]))]
    }
    data$Outcome[which(eval(com[i]))]=cf.data[,paste("Y",Valuev[i],sep=".")][which(eval(com[i]))]
  }

  return(list(data=data,cf.data=cf.data,coeff=coeff,Sigma=Sigma))
}

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
