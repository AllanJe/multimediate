#'
#'
#'@import stats
#'
#'
#'

CorCond = function(e,lmodel.m){
  NM=length(lmodel.m)
  k=1

  out=list()
  while(k < NM){
    for(l in (k+1):NM){
      a=paste("cor",names(lmodel.m[[k]]$model)[1],names(lmodel.m[[l]]$model)[1],sep=".")
      #print(a)
      out[[a]]=CorCond2(e,lmodel=list(lmodel.m[[k]],lmodel.m[[l]]))
    }
    k=k+1
  }
  Variance=rep(1,NM)
  for (z in 1:NM){
    if (inherits(lmodel.m[[z]],"lm") & !inherits(lmodel.m[[z]],"glm")){
      Variance[z]=var(lmodel.m[[z]]$residuals)
    }
  }
  sigmaestim=diag(Variance)
  correstim=diag(1,NM)
  i=1
  l=1
  vcor=c()
  vcov=c()
  while (i <NM){
    j=i+1
    for(k in j:NM){
      sigmaestim[i,k]=sigmaestim[k,i]=weighted.mean(out[[l]]$Covariance.estim,out[[l]]$Proportion,na.rm = TRUE)
      vcov=c(vcov,weighted.mean(out[[l]]$Covariance.estim,out[[l]]$Proportion,na.rm = TRUE))
      correstim[i,k]=correstim[k,i]=weighted.mean(out[[l]]$Correlation.estim,out[[l]]$Proportion,na.rm = TRUE)
      vcor=c(vcor,weighted.mean(out[[l]]$Correlation.estim,out[[l]]$Proportion,na.rm = TRUE))
      l=l+1
    }
    i=i+1
  }
  out[["sigmaestim"]]=sigmaestim
  out[["vcov"]]=vcov
  out[["correstim"]]=correstim
  out[["vcor"]]=vcor
  return(out)
}


#' @importFrom rmutil int

CorCond2 = function(e,lmodel){

  is.Lm.Mi=inherits(lmodel[[1]],"lm") & !inherits(lmodel[[1]],"glm")
  is.Lm.Mj=inherits(lmodel[[2]],"lm") & !inherits(lmodel[[2]],"glm")


  if(is.Lm.Mi & is.Lm.Mj){
    Estim=data.frame(NA)
    Estim$Proportion=1
    Estim$Covariance.estim=cov(lmodel[[1]]$residuals,lmodel[[2]]$residuals)
    Estim$Correlation.estim=cor(lmodel[[1]]$residuals,lmodel[[2]]$residuals)
  }

  else{
    Mi=names(lmodel[[1]]$model)[1]
    Mj=names(lmodel[[2]]$model)[1]
    data=cbind(lmodel[[1]]$model,lmodel[[2]]$model)[,unique(c(names(lmodel[[1]]$model),names(lmodel[[2]]$model)))]
    datav=data[,-c(which(names(data)==Mi),which(names(data)==Mj))]
    value=unique(datav)

    if(is.null(dim(value))){
      value=data.frame(t(t(value)))
      colnames(value)=names(lmodel[[1]]$model)[2]
      datav=data.frame(t(t(datav)))
      colnames(datav)=names(lmodel[[1]]$model)[2]
    }
    else{
      for (y in 1:dim(value)[2]){

        if (is.factor(value[,y])){
          test=length(grep(pattern = "[a-z][1-9][a-z]", value[,y], value = TRUE, fixed=FALSE))
          if (test==0){
            value[,y]=as.numeric(as.character(value[,y]))
          }

        }
      }
    }

    is.Polr.Mi=inherits(lmodel[[1]],"polr")
    is.Polr.Mj=inherits(lmodel[[2]],"polr")
    is.Glm.Mi=inherits(lmodel[[1]],"glm")
    is.Glm.Mj=inherits(lmodel[[2]],"glm")

    if(is.Polr.Mi & is.Polr.Mj){
      p=sort(as.numeric(as.character(unique(data[,Mi]))))
      interceptp=as.vector(c(-Inf,lmodel[[1]]$zeta,Inf))
      coefp=as.vector(lmodel[[1]]$coefficients)
      q=sort(as.numeric(as.character(unique(data[,Mj]))))
      interceptq=as.vector(c(-Inf,lmodel[[2]]$zeta,Inf))
      coefq=as.vector(lmodel[[2]]$coefficients)
    }

    if(is.Glm.Mi & is.Glm.Mj){
      p=c(0,1)
      interceptp=c(-Inf,-summary(lmodel[[1]])$coefficients[1,1],Inf)
      coefp=as.vector(summary(lmodel[[1]])$coefficients[-1,1])
      q=c(0,1)
      interceptq=c(-Inf,-summary(lmodel[[2]])$coefficients[1,1],Inf)
      coefq=as.vector(summary(lmodel[[2]])$coefficients[-1,1])
    }

    else if(is.Polr.Mi & is.Glm.Mj){
      p=sort(as.numeric(as.character(unique(data[,Mi]))))
      interceptp=as.vector(c(-Inf,lmodel[[1]]$zeta,Inf))
      coefp=as.vector(lmodel[[1]]$coefficients)
      q=c(0,1)
      interceptq=c(-Inf,-summary(lmodel[[2]])$coefficients[1,1],Inf)
      coefq=as.vector(summary(lmodel[[2]])$coefficients[-1,1])
    }

    else if(is.Glm.Mi & is.Polr.Mj){
      p=c(0,1)
      interceptp=c(-Inf,-summary(lmodel[[1]])$coefficients[1,1],Inf)
      coefp=as.vector(summary(lmodel[[1]])$coefficients[-1,1])
      q=sort(as.numeric(as.character(unique(data[,Mj]))))
      interceptq=as.vector(c(-Inf,lmodel[[2]]$zeta,Inf))
      coefq=as.vector(lmodel[[2]]$coefficients)
    }

    else if(is.Polr.Mi & is.Lm.Mj){
      p=sort(as.numeric(as.character(unique(data[,Mi]))))
      interceptp=as.vector(c(-Inf,lmodel[[1]]$zeta,Inf))
      coefp=as.vector(lmodel[[1]]$coefficients)
    }

    else if(is.Lm.Mi & is.Polr.Mj){
      p=sort(as.numeric(as.character(unique(data[,Mj]))))
      interceptp=as.vector(c(-Inf,lmodel[[2]]$zeta,Inf))
      coefp=as.vector(lmodel[[2]]$coefficients)
    }

    else if(is.Glm.Mi & is.Lm.Mj){
      p=c(0,1)
      interceptp=c(-Inf,-summary(lmodel[[1]])$coefficients[1,1],Inf)
      coefp=as.vector(summary(lmodel[[1]])$coefficients[-1,1])
    }

    else if (is.Lm.Mi & is.Glm.Mj){
      p=c(0,1)
      interceptp=c(-Inf,-summary(lmodel[[2]])$coefficients[1,1],Inf)
      coefp=as.vector(summary(lmodel[[2]])$coefficients[-1,1])
    }

    Value=NULL
    Proportion=NULL
    Covariance.estim=NULL
    Correlation.estim=NULL
    for (i in 1:dim(value)[1]){
      com = parse(text= paste(paste(names(datav),value[i,], sep = "=="), collapse = " & "))
      datai=subset(data,eval(com))
      valueMi=c(1,as.numeric(as.character(value[i,names(value) %in% names(lmodel[[1]]$model)])))
      valueMj=c(1,as.numeric(as.character(value[i,names(value) %in% names(lmodel[[2]]$model)])))

      covMiMj=cov(as.numeric(as.character(datai[,Mi])),as.numeric(as.character(datai[,Mj])))

      if( (is.Polr.Mi & is.Polr.Mj) | (is.Glm.Mi & is.Polr.Mj) | (is.Polr.Mi & is.Glm.Mj) | (is.Glm.Mi & is.Glm.Mj) ){
        bornep=borneq=NULL
        for (o in 1:length(interceptp)){
          bornep=c(bornep,sum(c(interceptp[o],-coefp)*valueMi))
        }
        for (o in 1:length(interceptq)){
          borneq=c(borneq,sum(c(interceptq[o],-coefq)*valueMj))
        }

        pMit=pnorm(bornep)
        pMit=pMit[2:length(pMit)]-pMit[1:(length(pMit)-1)]
        pMi=sum(p*pMit)
        pMjt=pnorm(borneq)
        pMjt=pMjt[2:length(pMjt)]-pMjt[1:(length(pMjt)-1)]
        pMj=sum(q*pMjt)
        cor=DichotomiePD(fPD,e,covMiMj,pMi,bornep,p,pMj,borneq,q)
      }

      if ((is.Polr.Mi & is.Lm.Mj) | (is.Lm.Mi & is.Polr.Mj) | (is.Glm.Mi & is.Lm.Mj) | (is.Lm.Mi & is.Glm.Mj))  {
        if(is.Polr.Mi & is.Lm.Mj){
          valueM=valueMi
        }

        else if(is.Lm.Mi & is.Polr.Mj){
          valueM=valueMj
        }

        else if(is.Glm.Mi & is.Lm.Mj){
          valueM=valueMi
        }

        else if (is.Lm.Mi & is.Glm.Mj){
          valueM=valueMj
        }
        bornep=NULL
        for (o in 1:length(interceptp)){
          bornep=c(bornep,sum(c(interceptp[o],-coefp)*valueM))
        }
        g=function(x){
          #out=sqrt((1-rho^2)/(1+rho^2))*x*dnorm(x,0,(1+rho^2)/(1-rho^2))
          return(x*dnorm(x,0,1))
        }
        smp=NULL
        for (e in 1:length(p)){
          smp=c(smp,p[e]*int(g, a=bornep[e], b=bornep[e+1]))
        }
        cor=covMiMj/sum(smp) #DichotomieCP(fCP,e,covMiMj,bornep,p)
      }


      Valuetf1=paste(value[i,])
      Valuetf2=NULL
      for (j in 1:length(Valuetf1)){
        Valuetf2=paste(Valuetf2,Valuetf1[j],sep="")
      }
      Value=c(Value,Valuetf2)
      Proportion=c(Proportion,dim(datai)[1]/dim(data)[1])
      Covariance.estim=c(Covariance.estim,round(covMiMj,4))
      Correlation.estim=c(Correlation.estim,round(cor,4))
    }

    Estim=data.frame(Value,Proportion,Covariance.estim,Correlation.estim)

  }

  return(Estim)
}

DichotomiePD = function(fPD,e,covMiMj,pMi,bornep,p,pMj,borneq,q){
  u=-1
  v=1
  K=0
  while (abs(u-v)>=e/2) {
    m=(u+v)/2
    fu=fPD(rho=u,covMiMj,pMi,bornep,p,pMj,borneq,q)
    fm=fPD(rho=m,covMiMj,pMi,bornep,p,pMj,borneq,q)
    if (fu*fm<0)
    {v=m}
    else
    {u=m}
    K=K+1
  }
  return(u)
}

#' @importFrom mvtnorm pmvnorm
fPD = function(rho,covMiMj,pMi,bornep,p,pMj,borneq,q){
  Sigma1=matrix(c(1,rho,rho,1),2,2)
  pmvn=NULL
  for (i in 1:length(p)){
    for(j in 1:length(q)){
      pmvn=c(pmvn,p[i]*q[j]*pmvnorm(lower=c(bornep[i],borneq[j]),upper=c(bornep[i+1],borneq[j+1]),sigma=Sigma1))
    }
  }
  out=covMiMj-sum(pmvn)+pMi*pMj
  return(out)
}
