#' multimediate
#'
#'
#' 'multimediate' is used to estimate various quantities for causal mediation analysis, including average causal mediation effects (indirect effect), average direct effects, proportions mediated, and total effect, in presence of multiple mediators uncausally related. With a binary outcome, 'multimediate' also estimate average causal mediation effects on OR scale and logOR scale.
#'
#'@param lmodel.m list of fitted models object for mediators. Can be of class 'lm', 'polr','glm'.
#'@param correlated a logical value. if 'FALSE' a identity matrix is used for the matrix of correlation of mediators; if 'TRUE' matrix of correlation is estimated. Default is 'FALSE'.
#'@param model.y a fitted model object for the outcome. Can be of class 'lm', 'polr','glm'.
#'@param treat a character string indicating the name of the treatment variable used in the models. The treatment can be either binary (integer or a two-valued factor) or continuous (numeric).
#'@param treat.value value of the treatment variable used as the treatment condition. Default is 1.
#'@param control.value value of the treatment variable used as the control condition. Default is 0.
#'@param J number of Monte Carlo draws for quasi-Bayesian approximation.
#'@param conf.level level of the returned two-sided confidence intervals. Default is to return the 2.5 and 97.5 percentiles of the simulated quantities.
#'@param fun the function used to compute the point estimate of the effects of interest from its empirical distribution. The function mean or median can be used. Default is the function mean.
#'
#'@return multimediate returns an object of class "mm", a list that contains the components listed below.
#' The function summary (i.e., summary.mm) can be used to obtain a table of the results.
#' @export
#'
#' @importFrom MASS mvrnorm
#' @importFrom mvtnorm rmvnorm
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar

multimediate=function(lmodel.m,correlated=FALSE,model.y,treat,treat.value=1,control.value=0,J=1000,conf.level=0.95,fun=mean){

  N=dim(lmodel.m[[1]]$model)[1]
  NM=length(lmodel.m)

  for(nm in 1:NM){
    testmm=inherits(lmodel.m[[nm]], c("lm","glm","polr"))
    if(testmm==FALSE){
      stop("Mediator's model not implemented yet.")
    }

    if (!inherits(lmodel.m[[nm]],"lm")){
      for(c in 1:dim(lmodel.m[[nm]]$model)[2]){
        if(!is.factor(lmodel.m[[nm]]$model[,c])){
          stop("covariate ",names(lmodel.m[[nm]]$model)[c]," in mediators models must be discretize and/or factorize.")
        }
      }


    }
  }

  testy=inherits(model.y, c("lm","glm","polr"))
  if(testy==FALSE){
    stop("Outcome's model not implemented yet.")
  }

  mediator=c()
  MModel = list()
  for (nm in 1:NM){
    MModel[[nm]]=rmvnorm(J, mean = c(coef(lmodel.m[[nm]]),lmodel.m[[nm]]$zeta), sigma = vcov(lmodel.m[[nm]]))
    #MModel[[nm]]=rmvnorm(J, mean = c(coef(lmodel.m[[nm]])), sigma = vcov(lmodel.m[[nm]]))
    mediator=c(mediator,names(lmodel.m[[nm]]$model)[1])
  }

  if (correlated==TRUE){
    mcov=CorCond(e=10^(-10),lmodel.m)
    error=mvrnorm(n=N*J,mu=rep(0,NM),mcov$sigmaestim,tol=100)
    cov=mcov$sigmaestim
    cor=mcov$correstim
  }
  else{
    error=mvrnorm(n=N*J,mu=rep(0,NM),diag(rep(1,NM)))
    cov=diag(rep(1,NM))
  }
  YModel = rmvnorm(J, mean = c(coef(model.y),model.y$zeta), sigma = vcov(model.y))


  if ( !is.null(model.y$family$link)){
    if (model.y$family$link=="logit"){
      Yerror <- rlogis(N,0,1)
    }
    else{
      Yerror <- rnorm(N,0,1)
    }
  }
  else {
    Yerror <- rnorm(N,0,1)
  }
  PredictM1<-PredictM0<-PredictM1b<-PredictM0b<- array(0, dim=c(J,N,NM))


  print("Simulation of counterfactuals mediators")
  pb <- txtProgressBar(min = 0, max = NM, style = 3,title ="Simulation of counterfactuals mediators")
  for (nm in 1:NM){
    pred.data.t <- pred.data.c <- model.frame(lmodel.m[[nm]])

    if (is.factor(lmodel.m[[nm]]$model[,treat])) {
      pred.data.t[, treat] <- factor(treat.value, levels = levels(lmodel.m[[nm]]$model[, treat]))
      pred.data.c[, treat] <- factor(control.value,levels= levels(lmodel.m[[nm]]$model[, treat]))
    }
    else {
      pred.data.t[, treat] <- treat.value
      pred.data.c[, treat] <- control.value
    }


    if(inherits(lmodel.m[[nm]],"polr")){
      mmat.t <- model.matrix(terms(lmodel.m[[nm]]), data = pred.data.t)
      mmat.c <- model.matrix(terms(lmodel.m[[nm]]), data = pred.data.c)
      mmat.t=mmat.t[,-1]
      mmat.c=mmat.c[,-1]

      if (is.null(dim(mmat.t))){
        muM1 <- tcrossprod(MModel[[nm]][,1], mmat.t)
        muM0 <- tcrossprod(MModel[[nm]][,1], mmat.c)
      }
      else{
        muM1 <- tcrossprod(MModel[[nm]][,1:dim(mmat.t)[2]], mmat.t)
        muM0 <- tcrossprod(MModel[[nm]][,1:dim(mmat.t)[2]], mmat.c)
      }


      PredictM1[,,nm] <- PredictM1b[,,nm] <- array(muM1,dim=c(J,N)) + array(error[,nm], dim=c(J,N))
      PredictM0[,,nm] <- PredictM0b[,,nm] <- array(muM0,dim=c(J,N)) + array(error[,nm], dim=c(J,N))

      seuil=c(-Inf,lmodel.m[[nm]]$zeta,Inf)
      # if (is.null(dim(mmat.t))){
      #   seuil=cbind(-Inf,MModel[[nm]][,-1],Inf)}
      # else{
      #   seuil=cbind(-Inf,MModel[[nm]][,-(1:dim(mmat.t)[2])],Inf)
      # }
      for (k in 1:length(lmodel.m[[nm]]$lev)){
        for (j in 1:J){
          #for (n in 1:N){
          # a=which(PredictM1b[,n,nm]>seuil[,k] & PredictM1b[,n,nm]<=seuil[,k+1])
          # b=which(PredictM0b[,n,nm]>seuil[,k] & PredictM0b[,n,nm]<=seuil[,k+1])
            # a=which(PredictM1b[,n,nm]>seuil[k] & PredictM1b[,n,nm]<=seuil[k+1])
            # b=which(PredictM0b[,n,nm]>seuil[k] & PredictM0b[,n,nm]<=seuil[k+1])
          a=which(PredictM1b[j,,nm]>seuil[k] & PredictM1b[j,,nm]<=seuil[k+1])
          b=which(PredictM0b[j,,nm]>seuil[k] & PredictM0b[j,,nm]<=seuil[k+1])
          # PredictM1[a,n,nm]=lmodel.m[[nm]]$lev[k]
          # PredictM0[b,n,nm]=lmodel.m[[nm]]$lev[k]
          PredictM1[j,a,nm]=lmodel.m[[nm]]$lev[k]
          PredictM0[j,b,nm]=lmodel.m[[nm]]$lev[k]
        }
      }
    }
    else{
      mmat.t <- model.matrix(terms(lmodel.m[[nm]]), data = pred.data.t)
      mmat.c <- model.matrix(terms(lmodel.m[[nm]]), data = pred.data.c)

      muM1 <- tcrossprod(MModel[[nm]], mmat.t)
      muM0 <- tcrossprod(MModel[[nm]], mmat.c)

      PredictM1[,,nm] <- array(muM1,dim=c(J,N)) + array(error[,nm], dim=c(J,N))
      PredictM0[,,nm] <- array(muM0,dim=c(J,N)) + array(error[,nm], dim=c(J,N))
      if (inherits(lmodel.m[[nm]],"glm")){
        PredictM1[,,nm]=(PredictM1[,,nm]>0)*1
        PredictM0[,,nm]=(PredictM0[,,nm]>0)*1
      }
    }

    setTxtProgressBar(pb, nm,title)
  }
  close(pb)


  effect.tmp.NM=OR.NM=array(NA, dim = c(N, J, 2, NM))
  effect.tmp=OR=array(NA, dim = c(N, J, 4))
  #OR.NM=array(NA, dim = c(J, 2, NM))
  #OR=array(NA, dim = c(J, 4))
  #OR.NM.polr=array(NA, dim = c(J, 2, NM,(length(model.y$lev)-2)))
  #OR.polr=array(NA, dim = c(J, 4,(length(model.y$lev)-2)))

  title=paste("Simulation of counterfactuals outcomes ",1:4,"/4",sep="")
  for (e in 1:4) {


    tt <- switch(e, c(1, 1, 1, 0), c(0, 0, 1, 0),
                 c(1, 0, 1, 1), c(1, 0, 0, 0))
    Pr0<-Pr1<-ORPr0<-ORPr1<- matrix(NA, nrow = N, ncol = J)
    if (NM!=1 & e<=2) {
      Pr0.NM<-Pr1.NM <-ORPr0.NM<-ORPr1.NM <- array(NA,dim=c(N,J,NM))
    }
    print(title[e])
    pb <- txtProgressBar(min = 0, max = J, style = 3,title=title[e])
    for (j in 1:J) {
      pred.data.t <- pred.data.c <-model.frame(model.y)

      cat.t <- ifelse(tt[1], treat.value, control.value)
      cat.c <- ifelse(tt[2], treat.value, control.value)

      if (is.factor(model.y$model[, treat])) {
        pred.data.t[, treat] <- factor(cat.t, levels = levels(model.y$model[, treat]))
        pred.data.c[, treat] <- factor(cat.c, levels = levels(model.y$model[, treat]))

      }
      else {
        pred.data.t[, treat] <- cat.t
        pred.data.c[, treat] <- cat.c

      }


      if (tt[3]==1){
        PredictMt=PredictM1[j,,]
      }
      else {
        PredictMt=PredictM0[j,,]
      }
      if (tt[4]==1){
        PredictMc=PredictM1[j,,]
      }
      else {
        PredictMc=PredictM0[j,,]
      }

      pred.data.t[, mediator] <- PredictMt
      pred.data.c[, mediator] <- PredictMc




      for (nm in 1:NM){
        if(inherits(lmodel.m[[nm]],"polr")){
          pred.data.c[,mediator[nm]]=as.factor(pred.data.c[,mediator[nm]])
          levels(pred.data.c[,mediator[nm]])=lmodel.m[[nm]]$lev
          pred.data.t[,mediator[nm]]=as.factor(pred.data.t[,mediator[nm]])
          levels(pred.data.t[,mediator[nm]])=lmodel.m[[nm]]$lev
        }

        if(inherits(lmodel.m[[nm]],"lm") & !inherits(lmodel.m[[nm]],"glm")){
          pred.data.c[,mediator[nm]]=as.numeric(pred.data.c[,mediator[nm]])
          pred.data.t[,mediator[nm]]=as.numeric(pred.data.t[,mediator[nm]])
        }

        if(inherits(lmodel.m[[nm]],"glm")){
          pred.data.c[,mediator[nm]]=as.factor(pred.data.c[,mediator[nm]])
          levels(pred.data.c[,mediator[nm]])=levels(lmodel.m[[nm]]$model[,1])
          pred.data.t[,mediator[nm]]=as.factor(pred.data.t[,mediator[nm]])
          levels(pred.data.t[,mediator[nm]])=levels(lmodel.m[[nm]]$model[,1])
        }
      }

      ymat.t=ymat.c=model.matrix(model.y)
      trans.t <- model.matrix(terms(model.y), data = pred.data.t)
      ymat.t[,colnames(ymat.t) %in% colnames(trans.t)]=trans.t
      ymat.t[,!(colnames(ymat.t) %in% colnames(trans.t))]=0

      trans.c <- model.matrix(terms(model.y), data = pred.data.c)
      ymat.c[,colnames(ymat.c) %in% colnames(trans.c)]=trans.c
      ymat.c[,!(colnames(ymat.c) %in% colnames(trans.c))]=0

      if (!inherits(model.y,"polr")){
        Pr1[,j] <- t(as.matrix(YModel[j, ])) %*% t(ymat.t) + Yerror
        Pr0[,j] <- t(as.matrix(YModel[j, ])) %*% t(ymat.c) + Yerror
        ORPr1[,j] <- t(as.matrix(YModel[j, ])) %*% t(ymat.t)
        ORPr0[,j] <- t(as.matrix(YModel[j, ])) %*% t(ymat.c)
      }
      else{
        ymat.t=ymat.t[,-1]
        ymat.c=ymat.c[,-1]
        Pr1[,j] <- t(as.matrix(YModel[j,1:dim(ymat.t)[2]])) %*% t(ymat.t)+ Yerror
        Pr0[,j] <- t(as.matrix(YModel[j,1:dim(ymat.c)[2]])) %*% t(ymat.c)+ Yerror

        #seuilY=cbind(-Inf,YModel[,-(1:dim(ymat.t)[2]+1)],Inf)
        seuilY=c(-Inf,model.y$zeta,Inf)
        Pr1b=as.numeric(Pr1[,j])
        Pr0b=as.numeric(Pr0[,j])
        for (k in 1:length(model.y$lev)){
          #a=which(Pr1b>seuilY[j,k] & Pr1b<=seuilY[j,k+1])
          #b=which(Pr0b>seuilY[j,k] & Pr0b<=seuilY[j,k+1])
          a=which(Pr1b>seuilY[k] & Pr1b<=seuilY[k+1])
          b=which(Pr0b>seuilY[k] & Pr0b<=seuilY[k+1])
          Pr1[a,j]=model.y$lev[k]
          Pr0[b,j]=model.y$lev[k]
        }
        Pr1=matrix(as.numeric(Pr1),N,J)
        Pr0=matrix(as.numeric(Pr0),N,J)
      }

      if (NM!=1 & e<=2) {

        for (nm in 1:NM){
          if (tt[3]==1){
            PredictMt=PredictM1[j,,nm]
          }
          else {
            PredictMt=PredictM0[j,,nm]
          }
          if (tt[1]==1){
            PredictWt=PredictM1[j,,-nm]
          }
          else {
            PredictWt=PredictM0[j,,-nm]
          }

          if (tt[4]==1){
            PredictMc=PredictM1[j,,nm]
          }
          else {
            PredictMc=PredictM0[j,,nm]
          }
          if (tt[2]==1){
            PredictWc=PredictM1[j,,-nm]
          }
          else {
            PredictWc=PredictM0[j,,-nm]
          }


          pred.data.t[, mediator[nm]] <- PredictMt
          pred.data.t[, mediator[-nm]] <- PredictWt
          pred.data.c[, mediator[nm]] <- PredictMc
          pred.data.c[, mediator[-nm]] <- PredictWc
          for (nmt in 1:NM){
            if(inherits(lmodel.m[[nmt]],"polr")){

              pred.data.c[,mediator[nmt]]=as.factor(pred.data.c[,mediator[nmt]])
              levels(pred.data.c[,mediator[nmt]])=lmodel.m[[nmt]]$lev
              pred.data.t[,mediator[nmt]]=as.factor(pred.data.t[,mediator[nmt]])
              levels(pred.data.t[,mediator[nmt]])=lmodel.m[[nmt]]$lev
            }

            if(inherits(lmodel.m[[nmt]],"lm") & !inherits(lmodel.m[[nmt]],"glm")){
              pred.data.c[,mediator[nmt]]=as.numeric(pred.data.c[,mediator[nmt]])
              pred.data.t[,mediator[nmt]]=as.numeric(pred.data.t[,mediator[nmt]])
            }

            if(inherits(lmodel.m[[nmt]],"glm")){
              pred.data.c[,mediator[nmt]]=as.factor(pred.data.c[,mediator[nmt]])
              levels(pred.data.c[,mediator[nmt]])=levels(lmodel.m[[nmt]]$model[,1])
              pred.data.t[,mediator[nmt]]=as.factor(pred.data.t[,mediator[nmt]])
              levels(pred.data.t[,mediator[nmt]])=levels(lmodel.m[[nmt]]$model[,1])
            }
          }

          ymat.t=ymat.c=model.matrix(model.y)
          trans.t <- model.matrix(terms(model.y), data = pred.data.t)
          ymat.t[,colnames(ymat.t) %in% colnames(trans.t)]=trans.t
          ymat.t[,!(colnames(ymat.t) %in% colnames(trans.t))]=0

          trans.c <- model.matrix(terms(model.y), data = pred.data.c)
          ymat.c[,colnames(ymat.c) %in% colnames(trans.c)]=trans.c
          ymat.c[,!(colnames(ymat.c) %in% colnames(trans.c))]=0


          if(!inherits(model.y,"polr")){
            Pr1.NM[,j,nm] <- t(as.matrix(YModel[j, ])) %*% t(ymat.t)+ Yerror
            Pr0.NM[,j,nm] <- t(as.matrix(YModel[j, ])) %*% t(ymat.c)+ Yerror
            ORPr1.NM[,j,nm] <- t(as.matrix(YModel[j, ])) %*% t(ymat.t)
            ORPr0.NM[,j,nm] <- t(as.matrix(YModel[j, ])) %*% t(ymat.c)
          }
          else{
            ymat.t=ymat.t[,-1]
            ymat.c=ymat.c[,-1]
            Pr1.NM[,j,nm] <- t(as.matrix(YModel[j,1:dim(ymat.t)[2]])) %*% t(ymat.t)+ Yerror
            Pr0.NM[,j,nm] <- t(as.matrix(YModel[j,1:dim(ymat.c)[2]])) %*% t(ymat.c)+ Yerror

            Pr1b.NM=as.numeric(Pr1.NM[,j,nm])
            Pr0b.NM=as.numeric(Pr0.NM[,j,nm])
            for (k in 1:length(model.y$lev)){
              # a=which(Pr1b.NM>seuilY[j,k] & Pr1b.NM<=seuilY[j,k+1])
              # b=which(Pr0b.NM>seuilY[j,k] & Pr0b.NM<=seuilY[j,k+1])
              a=which(Pr1b.NM>seuilY[k] & Pr1b.NM<=seuilY[k+1])
              b=which(Pr0b.NM>seuilY[k] & Pr0b.NM<=seuilY[k+1])
              Pr1.NM[a,j,nm]=model.y$lev[k]
              Pr0.NM[b,j,nm]=model.y$lev[k]
            }
            Pr1.NM=array(as.numeric(Pr1.NM),dim=c(N,J,NM))
            Pr0.NM=array(as.numeric(Pr0.NM),dim=c(N,J,NM))
          }
        }
      }
      setTxtProgressBar(pb, j,title=title[e])
    }
    close(pb)

    if (!is.null(model.y$family)){

      effect.tmp[,,e]=(Pr1>0)-(Pr0>0)
      if (model.y$family$link=="logit"){

        expit=function(x){
          res=1/(1+exp(-x))
          return(res)
        }

        OR[,,e]=((1-expit(ORPr0))/expit(ORPr0))*((expit(ORPr1))/(1-expit(ORPr1)))
        #OR[,e]=(1-apply(Pr0>0,c(2),mean))/(1-apply(Pr1>0,c(2),mean))*(1+(apply(effect.tmp[,,e],c(2),mean)/apply(Pr0>0,c(2),mean)))
      }
      if (NM!=1 & e<=2){

        effect.tmp.NM[,,e,]=(Pr1.NM>0)-(Pr0.NM>0)
        if (model.y$family$link=="logit"){
          OR.NM[,,e,]=((1-expit(ORPr0.NM))/expit(ORPr0.NM))*((expit(ORPr1.NM))/(1-expit(ORPr1.NM)))
          #OR.NM[,e,]=(1-apply(Pr0.NM>=0, c(2,3), mean))/(1-apply(Pr1.NM>=0, c(2,3), mean))*(1+(apply(effect.tmp.NM[,,e,],c(2,3),mean)/apply(Pr0.NM>=0, c(2,3), mean)))
        }

      }
    }

    else {
      effect.tmp[,,e]=Pr1-Pr0
      # if (inherits(model.y,"polr")){
      #   for( pol in 2:(length(model.y$lev)-1)){
      #     print(pol)
      #     lev=as.numeric(model.y$lev)[pol]
      #     OR.polr[,e,pol]=(apply(Pr0<=lev,c(2),mean))/(apply(Pr1<=lev,c(2),mean))*(apply(Pr1>lev,c(2),mean))/(apply(Pr0>lev,c(2),mean))
      #   }
      #
      # }
      if (NM!=1 & e<=2){
        effect.tmp.NM[,,e,]=Pr1.NM-Pr0.NM
        #  if (inherits(model.y,"polr")){
        #   for( pol in 2:(length(model.y$lev)-1)){
        #     print(pol)
        #     lev=as.numeric(model.y$lev)[pol]
        #     OR.NM.polr[,e,,pol]=(apply(Pr0.NM<=lev,c(2,3),mean))/(apply(Pr1.NM<=lev,c(2,3),mean))*(apply(Pr1.NM>lev,c(2,3),mean))/(apply(Pr0.NM>lev,c(2,3),mean))
        #   }
        # }
      }
    }
}
  print("Computing average point estimates together with p-values and confidence intervals")
  # Step 3.3 : Compute the effects.
  # Compute the average effects. That is, we simply take the difference accross
  # two outcome predictions under treatment and the two outcome predictions under
  # control and average across the predictions for each of the n units in the study.
  # This provides us with effect^(j)(t), which is one Monte Carlo draw of the average effect

  et1 <- effect.tmp[, , 1] # joint mediated effect delta(1)
  et2 <- effect.tmp[, , 2] # joint mediated effect delta(0)
  et3 <- effect.tmp[, , 3] # direct effect zeta(1)
  et4 <- effect.tmp[, , 4] # direct effect zeta(0)

  delta.1 <- t(as.matrix(apply(et1, 2, mean,na.rm=TRUE)))
  delta.0 <- t(as.matrix(apply(et2, 2, mean,na.rm=TRUE)))

  zeta.1 <- t(as.matrix(apply(et3, 2, mean,na.rm=TRUE)))
  zeta.0 <- t(as.matrix(apply(et4, 2, mean,na.rm=TRUE)))

  tau <- (zeta.1 + delta.0 + zeta.0 + delta.1)/2
  nu.0 <- delta.0/tau
  nu.1 <- delta.1/tau
  delta.avg <- (delta.1 + delta.0)/2
  zeta.avg <- (zeta.1 + zeta.0)/2
  nu.avg <- (nu.1 + nu.0)/2
  d0 <- fun(delta.0,na.rm=TRUE)
  d1 <- fun(delta.1,na.rm=TRUE)
  z1 <- fun(zeta.1,na.rm=TRUE)
  z0 <- fun(zeta.0,na.rm=TRUE)
  tau.coef <- fun(tau,na.rm=TRUE)
  n0 <- fun(nu.0*is.finite(nu.0),na.rm=TRUE)#median(nu.0,na.rm=TRUE) # when total effect is 0 it produce NaN
  n1 <- fun(nu.1*is.finite(nu.1),na.rm=TRUE)#median(nu.1,na.rm=TRUE)
  d.avg <- (d0 + d1)/2
  z.avg <- (z0 + z1)/2
  n.avg <- (n0 + n1)/2

  if (!is.null(model.y$family)){
    if ( model.y$family$link=="logit"){
      ORet1 <- OR[ , ,1] # joint mediated effect ORdelta(1)
      ORet2 <- OR[ , ,2] # joint mediated effect ORdelta(0)
      ORet3  <- OR[ , ,3] # direct effect ORzeta(1)
      ORet4  <- OR[ , ,4] # direct effect ORzeta(0)
      ORdelta.1 <- t(as.matrix(apply(ORet1, 2, mean,na.rm=TRUE)))
      ORdelta.0 <- t(as.matrix(apply(ORet2, 2, mean,na.rm=TRUE)))
      ORzeta.1 <- t(as.matrix(apply(ORet3, 2, mean,na.rm=TRUE)))
      ORzeta.0 <- t(as.matrix(apply(ORet4, 2, mean,na.rm=TRUE)))

      # ORdelta.1 <- OR[ , 1] # joint mediated effect ORdelta(1)
      # ORdelta.0 <- OR[ , 2] # joint mediated effect ORdelta(0)
      # ORzeta.1  <- OR[ , 3] # direct effect ORzeta(1)
      # ORzeta.0  <- OR[ , 4] # direct effect ORzeta(0)



      ORtau <- (ORzeta.1*ORdelta.0 + ORzeta.0*ORdelta.1)/2

      ORtau <- (ORzeta.1*ORdelta.0 + ORzeta.0*ORdelta.1)/2
      ORnu.1 <- ORzeta.0*(ORdelta.1-1)/(ORzeta.0*ORdelta.1-1)
      ORnu.0 <- ORzeta.1*(ORdelta.0-1)/(ORzeta.1*ORdelta.0-1)

      ORtau.coef <- fun(ORtau,na.rm=TRUE)#median(ORtau,na.rm=TRUE)

      ORdelta.avg <- (ORdelta.1 + ORdelta.0)/2
      ORnu.avg <- (ORnu.1 + ORnu.0)/2
      ORzeta.avg <- (ORzeta.1 + ORzeta.0)/2

      ORd1 <- fun(ORdelta.1*is.finite(ORdelta.1),na.rm=TRUE)#median(ORdelta.1*is.finite(ORdelta.1),na.rm=TRUE)
      ORn1 <- fun(ORnu.1*is.finite(ORnu.1),na.rm=TRUE)#median(ORnu.1*is.finite(ORnu.1),na.rm=TRUE)
      ORd0 <- fun(ORdelta.0*is.finite(ORdelta.0),na.rm=TRUE)#median(ORdelta.0*is.finite(ORdelta.0),na.rm=TRUE)
      ORn0 <- fun(ORnu.0*is.finite(ORnu.0),na.rm=TRUE)#median(ORnu.0*is.finite(ORnu.0),na.rm=TRUE)
      ORz1 <- fun(ORzeta.1*is.finite(ORzeta.1),na.rm=TRUE)#median(ORzeta.1*is.finite(ORzeta.1),na.rm=TRUE)
      ORz0 <- fun(ORzeta.0*is.finite(ORzeta.0),na.rm=TRUE)#median(ORzeta.0*is.finite(ORzeta.0),na.rm=TRUE)

      ORd.avg <- (ORd0 + ORd1)/2
      ORz.avg <- (ORz0 + ORz1)/2
      ORn.avg <- (ORn0 + ORn1)/2

      logORdelta.1 <- log(ORdelta.1)
      logORtau <- log(ORtau)
      logORnu.1 <- logORdelta.1/logORtau
      logORdelta.0 <- log(ORdelta.0)
      logORnu.0 <- logORdelta.0/logORtau
      logORzeta.1 <- log(ORzeta.1)
      logORzeta.0 <- log(ORzeta.0)

      logORtau.coef <- fun(logORtau*is.finite(logORtau),na.rm=TRUE)#median(logORtau*is.finite(logORtau),na.rm=TRUE)
      logORnu.avg <- (logORnu.1 + logORnu.0)/2
      logORn0 <- fun(logORnu.0*is.finite(logORnu.0),na.rm=TRUE)#median(logORnu.0*is.finite(logORnu.0),na.rm=TRUE)
      logORn1 <- fun(logORnu.1*is.finite(logORnu.1),na.rm=TRUE)#median(logORnu.1*is.finite(logORnu.1),na.rm=TRUE)
      logORn.avg <- (logORn0 + logORn1)/2

      logORdelta.avg <- (logORdelta.1 + logORdelta.0)/2
      logORzeta.avg <- (logORzeta.1 + logORzeta.0)/2
      logORd0 <- fun(logORdelta.0*is.finite(logORdelta.0),na.rm=TRUE)#median(logORdelta.0*is.finite(logORdelta.0),na.rm=TRUE)
      logORd1 <- fun(logORdelta.1*is.finite(logORdelta.1),na.rm=TRUE)#median(logORdelta.1*is.finite(logORdelta.1),na.rm=TRUE)
      logORz1 <- fun(logORzeta.1*is.finite(logORzeta.1),na.rm=TRUE)#median(logORzeta.1*is.finite(logORzeta.1),na.rm=TRUE)
      logORz0 <- fun(logORzeta.0*is.finite(logORzeta.0),na.rm=TRUE)#median(logORzeta.0*is.finite(logORzeta.0),na.rm=TRUE)
      logORtau.coef <- fun(logORtau*is.finite(logORtau),na.rm=TRUE)#median(logORtau*is.finite(logORtau),na.rm=TRUE)
      logORd.avg <- (logORd0 + logORd1)/2
      logORz.avg <- (logORz0 + logORz1)/2

    }}

  if (NM!=1){
    et5 <- effect.tmp.NM[, ,1,] # mediated effect delta(1)
    et6 <- effect.tmp.NM[, ,2,] # mediated effect delta(0)

    delta.1.NM <- t((apply(et5, c(2,3), mean)))
    eta.1.NM <- array(delta.1,dim=c(NM,J))-delta.1.NM
    delta.0.NM <- t(as.matrix(apply(et6, c(2,3), mean)))
    eta.0.NM <- array(delta.0,dim=c(NM,J))-delta.0.NM

    nu.0.NM <- delta.0.NM/array(tau,dim=c(NM,J))
    nu.1.NM <- delta.1.NM/array(tau,dim=c(NM,J))
    delta.avg.NM <- (delta.1.NM + delta.0.NM)/2
    nu.avg.NM <- (nu.1.NM + nu.0.NM)/2
    d0.NM <- apply(delta.0.NM*is.finite(delta.0.NM),1,fun,na.rm=TRUE)
    d1.NM <- apply(delta.1.NM*is.finite(delta.1.NM),1,fun,na.rm=TRUE)
    n0.NM <- apply(nu.0.NM*is.finite(nu.0.NM),1,fun,na.rm=TRUE)#apply(nu.0.NM,1,median,na.rm=TRUE)
    n1.NM <- apply(nu.1.NM*is.finite(nu.1.NM),1,fun,na.rm=TRUE)#apply(nu.1.NM,1,median,na.rm=TRUE)
    d.avg.NM <- (d0.NM + d1.NM)/2
    n.avg.NM <- (n0.NM + n1.NM)/2

    if (!is.null(model.y$family)){
      if (model.y$family$link=="logit"){
        ORet5 <- OR.NM[, ,1,] # mediated effect ORdelta(1)
        ORet6 <- OR.NM[, ,2,] # mediated effect ORdelta(0)
        ORdelta.1.NM <- t((apply(ORet5, c(2,3), mean)))
        OReta.1.NM <- array(ORdelta.1,dim=c(NM,J))/ORdelta.1.NM
        ORdelta.0.NM <- t(as.matrix(apply(ORet6, c(2,3), mean)))
        OReta.0.NM <- array(ORdelta.0,dim=c(NM,J))/ORdelta.0.NM

        ORnu.0.NM <- array(ORzeta.1,dim=c(NM,J))*OReta.0.NM*(ORdelta.0.NM-1)/(array(ORzeta.1*ORdelta.0,dim=c(NM,J))-1)
        ORnu.1.NM <- array(ORzeta.0,dim=c(NM,J))*OReta.1.NM*(ORdelta.1.NM-1)/(array(ORzeta.0*ORdelta.1,dim=c(NM,J))-1)
        ORdelta.avg.NM <- (ORdelta.1.NM + ORdelta.0.NM)/2
        ORnu.avg.NM <- (ORnu.1.NM + ORnu.0.NM)/2
        ORd0.NM <- apply(ORdelta.0.NM*is.finite(ORdelta.0.NM),1,fun,na.rm=TRUE)#apply(ORdelta.0.NM,2,median,na.rm=TRUE)
        ORd1.NM <- apply(ORdelta.1.NM*is.finite(ORdelta.1.NM),1,fun,na.rm=TRUE)#apply(ORdelta.1.NM,2,median,na.rm=TRUE)
        ORn0.NM <- apply(ORnu.0.NM*is.finite(ORnu.0.NM),1,fun,na.rm=TRUE)#apply(ORnu.0.NM,2,median,na.rm=TRUE)
        ORn1.NM <- apply(ORnu.1.NM*is.finite(ORnu.1.NM),1,fun,na.rm=TRUE)#apply(ORnu.1.NM,2,median,na.rm=TRUE)
        ORd.avg.NM <- (ORd0.NM + ORd1.NM)/2
        ORn.avg.NM <- (ORn0.NM + ORn1.NM)/2

        # ORnu.0.NM <- array(ORzeta.1,dim=c(J,NM))*(ORdelta.0.NM-1)/(array(ORzeta.1*ORdelta.0,dim=c(J,NM))-1)
        # ORnu.1.NM <- array(ORzeta.0,dim=c(J,NM))*(ORdelta.1.NM-1)/(array(ORzeta.0*ORdelta.1,dim=c(J,NM))-1)
        # ORdelta.avg.NM <- (ORdelta.1.NM + ORdelta.0.NM)/2
        # ORnu.avg.NM <- (ORnu.1.NM + ORnu.0.NM)/2
        # ORd0.NM <- apply(ORdelta.0.NM,2,mean,na.rm=TRUE)#apply(ORdelta.0.NM,2,median,na.rm=TRUE)
        # ORd1.NM <- apply(ORdelta.1.NM,2,mean,na.rm=TRUE)#apply(ORdelta.1.NM,2,median,na.rm=TRUE)
        # ORn0.NM <- apply(ORnu.0.NM,2,mean,na.rm=TRUE)#apply(ORnu.0.NM,2,median,na.rm=TRUE)
        # ORn1.NM <- apply(ORnu.1.NM,2,mean,na.rm=TRUE)#apply(ORnu.1.NM,2,median,na.rm=TRUE)
        # ORd.avg.NM <- (ORd0.NM + ORd1.NM)/2
        # ORn.avg.NM <- (ORn0.NM + ORn1.NM)/2

        logORdelta.1.NM <- log(ORdelta.1.NM)
        logORdelta.0.NM <- log(ORdelta.0.NM)
        logORnu.0.NM <- logORdelta.0.NM/array(logORtau,dim=c(NM,J))
        logORnu.1.NM <- logORdelta.1.NM/array(logORtau,dim=c(NM,J))

        logORnu.avg.NM <- (logORnu.1.NM + logORnu.0.NM)/2
        logORn0.NM <- apply(logORdelta.0.NM/logORtau.coef,1,fun,na.rm=TRUE)#apply(logORdelta.0.NM/logORtau.coef,2,median,na.rm=TRUE)
        logORn1.NM <- apply(logORdelta.1.NM/logORtau.coef,1,fun,na.rm=TRUE)#apply(logORdelta.1.NM/logORtau.coef,2,median,na.rm=TRUE)
        logORn.avg.NM <- (logORn0.NM + logORn1.NM)/2

        logORdelta.avg.NM <- (logORdelta.1.NM + logORdelta.0.NM)/2
        logORd0.NM <- apply(logORdelta.0.NM,1,fun,na.rm=TRUE)#apply(logORdelta.0.NM,2,median,na.rm=TRUE)
        logORd1.NM <- apply(logORdelta.1.NM,1,fun,na.rm=TRUE)#apply(logORdelta.1.NM,2,median,na.rm=TRUE)
        logORd.avg.NM <- (logORd0.NM + logORd1.NM)/2
      }}

  }
  else {
    delta.1.NM <- delta.0.NM <- NULL
    if (!is.null(model.y$family)){
      if (model.y$family$link=="logit"){
        ORdelta.1.NM <- ORdelta.0.NM <- logORdelta.1.NM <- logORdelta.0.NM <- NULL
      }}
  }


  low <- (1 - conf.level)/2
  high <- 1 - low

  d0.ci <- quantile(delta.0, c(low, high), na.rm = TRUE)
  d1.ci <- quantile(delta.1, c(low, high), na.rm = TRUE)
  tau.ci <- quantile(tau, c(low, high), na.rm = TRUE)
  z1.ci <- quantile(zeta.1, c(low, high), na.rm = TRUE)
  z0.ci <- quantile(zeta.0, c(low, high), na.rm = TRUE)
  n0.ci <- quantile(nu.0, c(low, high), na.rm = TRUE)
  n1.ci <- quantile(nu.1, c(low, high), na.rm = TRUE)
  d.avg.ci <- quantile(delta.avg, c(low, high), na.rm = TRUE)
  z.avg.ci <- quantile(zeta.avg, c(low, high), na.rm = TRUE)
  n.avg.ci <- quantile(nu.avg, c(low, high), na.rm = TRUE)

  d0.p <- pval(delta.0, d0)
  d1.p <- pval(delta.1, d1)
  d.avg.p <- pval(delta.avg, d.avg)
  z0.p <- pval(zeta.0, z0)
  z1.p <- pval(zeta.1, z1)
  z.avg.p <- pval(zeta.avg, z.avg)
  n0.p <- pval(nu.0, n0)
  n1.p <- pval(nu.1, n1)
  n.avg.p <- pval(nu.avg, n.avg)
  tau.p <- pval(tau, tau.coef)

  if (!is.null(model.y$family)){
    if (model.y$family$link=="logit"){
      ORd0.ci <- quantile(ORdelta.0, c(low, high), na.rm = TRUE)
      ORd1.ci <- quantile(ORdelta.1, c(low, high), na.rm = TRUE)
      ORtau.ci <- quantile(ORtau, c(low, high), na.rm = TRUE)
      ORz1.ci <- quantile(ORzeta.1, c(low, high), na.rm = TRUE)
      ORz0.ci <- quantile(ORzeta.0, c(low, high), na.rm = TRUE)
      ORn0.ci <- quantile(ORnu.0, c(low, high), na.rm = TRUE)
      ORn1.ci <- quantile(ORnu.1, c(low, high), na.rm = TRUE)

      ORd.avg.ci <- quantile(ORdelta.avg, c(low, high), na.rm = TRUE)
      ORz.avg.ci <- quantile(ORzeta.avg, c(low, high), na.rm = TRUE)
      ORn.avg.ci <- quantile(ORnu.avg, c(low, high), na.rm = TRUE)

      ORd0.p <- pval(ORdelta.0, ORd0,seu=1)
      ORd1.p <- pval(ORdelta.1, ORd1,seu=1)
      ORd.avg.p <- pval(ORdelta.avg, ORd.avg,seu=1)
      ORz0.p <- pval(ORzeta.0, ORz0,seu=1)
      ORz1.p <- pval(ORzeta.1, ORz1,seu=1)
      ORz.avg.p <- pval(ORzeta.avg, ORz.avg,seu=1)
      ORn0.p <- pval(ORnu.0, ORn0)
      ORn1.p <- pval(ORnu.1, ORn1)
      ORn.avg.p <- pval(ORnu.avg, ORn.avg)
      ORtau.p <- pval(ORtau, ORtau.coef,seu=1)

      logORn0.ci <- quantile(logORnu.0, c(low, high), na.rm = TRUE)
      logORn1.ci <- quantile(logORnu.1, c(low, high), na.rm = TRUE)
      logORn.avg.ci <- quantile(logORnu.avg, c(low, high), na.rm = TRUE)
      logORn0.p <- pval(logORnu.0, logORn0)
      logORn1.p <- pval(logORnu.1, logORn1)
      logORn.avg.p <- pval(logORnu.avg, logORn.avg)
      logORd0.ci <- quantile(logORdelta.0, c(low, high), na.rm = TRUE)
      logORd1.ci <- quantile(logORdelta.1, c(low, high), na.rm = TRUE)
      logORtau.ci <- quantile(logORtau, c(low, high), na.rm = TRUE)
      logORz1.ci <- quantile(logORzeta.1, c(low, high), na.rm = TRUE)
      logORz0.ci <- quantile(logORzeta.0, c(low, high), na.rm = TRUE)
      logORd.avg.ci <- quantile(logORdelta.avg, c(low, high), na.rm = TRUE)
      logORz.avg.ci <- quantile(logORzeta.avg, c(low, high), na.rm = TRUE)

      logORd0.p <- pval(logORdelta.0, logORd0)
      logORd1.p <- pval(logORdelta.1, logORd1)
      logORd.avg.p <- pval(logORdelta.avg, logORd.avg)
      logORz0.p <- pval(logORzeta.0, logORz0)
      logORz1.p <- pval(logORzeta.1, logORz1)
      logORz.avg.p <- pval(logORzeta.avg, logORz.avg)
      logORtau.p <- pval(logORtau, logORtau.coef)
    }}


  if(NM!=1){
    d0.ci.NM <- d1.ci.NM <- n0.ci.NM <- n1.ci.NM <- d.avg.ci.NM <- n.avg.ci.NM <- array(NA,dim=c(NM,2))
    d0.p.NM <- d1.p.NM <- d.avg.p.NM <- n0.p.NM <- n1.p.NM <- n.avg.p.NM <- rep(NA,NM)
    for (nm in 1:NM){
      d0.ci.NM[nm,] <- quantile(delta.0.NM[nm,], c(low, high), na.rm = TRUE)
      d1.ci.NM[nm,] <- quantile(delta.1.NM[nm,], c(low, high), na.rm = TRUE)
      n0.ci.NM[nm,] <- quantile(nu.0.NM[nm,], c(low, high), na.rm = TRUE)
      n1.ci.NM[nm,] <- quantile(nu.1.NM[nm,], c(low, high), na.rm = TRUE)
      d.avg.ci.NM[nm,] <- quantile(delta.avg.NM[nm,], c(low, high), na.rm = TRUE)
      n.avg.ci.NM[nm,] <- quantile(nu.avg.NM[nm,], c(low, high), na.rm = TRUE)

      d0.p.NM[nm] <- pval(delta.0.NM[nm,], d0.NM[nm])
      d1.p.NM[nm] <- pval(delta.1.NM[nm,], d1.NM[nm])
      d.avg.p.NM[nm] <- pval(delta.avg.NM[nm,], d.avg.NM[nm])
      n0.p.NM[nm] <- pval(nu.0.NM[nm,], n0.NM[nm])
      n1.p.NM[nm] <- pval(nu.1.NM[nm,], n1.NM[nm])
      n.avg.p.NM[nm] <- pval(nu.avg.NM[nm,], n.avg.NM[nm])
    }

    if (!is.null(model.y$family)){
      if (model.y$family$link=="logit"){
        ORd0.ci.NM <- ORd1.ci.NM <- ORn0.ci.NM <- ORn1.ci.NM <-ORd.avg.ci.NM <- ORn.avg.ci.NM <- array(NA,dim=c(NM,2))
        ORd0.p.NM <- ORd1.p.NM <- ORd.avg.p.NM <- ORn0.p.NM <- ORn1.p.NM <- ORn.avg.p.NM<- rep(NA,NM)

        logORd0.ci.NM <- logORd1.ci.NM <- logORn0.ci.NM <- logORn1.ci.NM <- logORd.avg.ci.NM <- logORn.avg.ci.NM <- array(NA,dim=c(NM,2))
        logORd0.p.NM <- logORd1.p.NM <- logORd.avg.p.NM <- logORn0.p.NM <- logORn1.p.NM <- logORn.avg.p.NM <- rep(NA,NM)
        # logORd0.p.NM <- logORd1.p.NM <- logORd.avg.p.NM <- logn0.p.NM <- logn1.p.NM <- logn.avg.p.NM <- rep(NA,NM)
        for (nm in 1:NM){
          ORd0.ci.NM[nm,] <- quantile(ORdelta.0.NM[nm,], c(low, high), na.rm = TRUE)
          ORd1.ci.NM[nm,] <- quantile(ORdelta.1.NM[nm,], c(low, high), na.rm = TRUE)
          ORn0.ci.NM[nm,] <- quantile(ORnu.0.NM[nm,], c(low, high), na.rm = TRUE)
          ORn1.ci.NM[nm,] <- quantile(ORnu.1.NM[nm,], c(low, high), na.rm = TRUE)
          ORd.avg.ci.NM[nm,] <- quantile(ORdelta.avg.NM[nm,], c(low, high), na.rm = TRUE)
          ORn.avg.ci.NM[nm,] <- quantile(ORnu.avg.NM[nm,], c(low, high), na.rm = TRUE)

          ORd0.p.NM[nm] <- pval(ORdelta.0.NM[nm,], ORd0.NM[nm],seu=1)
          ORd1.p.NM[nm] <- pval(ORdelta.1.NM[nm,], ORd1.NM[nm],seu=1)
          ORd.avg.p.NM[nm] <- pval(ORdelta.avg.NM[nm,], ORd.avg.NM[nm],seu=1)
          ORn0.p.NM[nm] <- pval(ORnu.0.NM[nm,], ORn0.NM[nm])
          ORn1.p.NM[nm] <- pval(ORnu.1.NM[nm,], ORn1.NM[nm])
          ORn.avg.p.NM[nm] <- pval(ORnu.avg.NM[nm,], ORn.avg.NM[nm])

          logORd0.ci.NM[nm,] <- quantile(logORdelta.0.NM[nm,], c(low, high), na.rm = TRUE)
          logORd1.ci.NM[nm,] <- quantile(logORdelta.1.NM[nm,], c(low, high), na.rm = TRUE)
          logORd.avg.ci.NM[nm,] <- quantile(logORdelta.avg.NM[nm,], c(low, high), na.rm = TRUE)
          logORn0.ci.NM[nm,] <- quantile(logORnu.0.NM[nm,], c(low, high), na.rm = TRUE)
          logORn1.ci.NM[nm,] <- quantile(logORnu.1.NM[nm,], c(low, high), na.rm = TRUE)
          logORn.avg.ci.NM[nm,] <- quantile(logORnu.avg.NM[nm,], c(low, high), na.rm = TRUE)
          logORd0.p.NM[nm] <- pval(logORdelta.0.NM[nm,], logORd0.NM[nm])
          logORd1.p.NM[nm] <- pval(logORdelta.1.NM[nm,], logORd1.NM[nm])
          logORd.avg.p.NM[nm] <- pval(logORdelta.avg.NM[nm,], logORd.avg.NM[nm])
          logORn0.p.NM[nm] <- pval(logORnu.0.NM[nm,], logORn0.NM[nm])
          logORn1.p.NM[nm] <- pval(logORnu.1.NM[nm,], logORn1.NM[nm])
          logORn.avg.p.NM[nm] <- pval(logORnu.avg.NM[nm,], logORn.avg.NM[nm])
        }


      }}
  }
  if (NM==1){
    out <- list(d0 = d0, d1 = d1, d0.ci = d0.ci, d1.ci = d1.ci, d0.p = d0.p, d1.p = d1.p, d0.sims = delta.0,d1.sims = delta.1,
                z0 = z0, z1 = z1, z0.ci = z0.ci,z1.ci = z1.ci, z0.p = z0.p, z1.p = z1.p, z0.sims = zeta.0, z1.sims = zeta.1,
                n0 = n0, n1 = n1, n0.ci = n0.ci, n1.ci = n1.ci, n0.p = n0.p, n1.p = n1.p, n0.sims = nu.0, n1.sims = nu.1,
                tau.coef = tau.coef, tau.ci = tau.ci, tau.p = tau.p, tau.sims = tau,
                d.avg = d.avg, d.avg.p = d.avg.p, d.avg.ci = d.avg.ci, d.avg.sims = delta.avg,
                z.avg = z.avg, z.avg.p = z.avg.p, z.avg.ci = z.avg.ci, z.avg.sims = zeta.avg,
                n.avg = n.avg, n.avg.p = n.avg.p, n.avg.ci = n.avg.ci, n.avg.sims = nu.avg,
                treat = treat, mediator = mediator,
                conf.level = conf.level,
                model.y = model.y, model.m = lmodel.m,
                control.value = control.value, treat.value = treat.value,
                nobs = N, sims = J,cov=cov,cor=cor)
    if (!is.null(model.y$family)){
      if (model.y$family$link=="logit"){
        out <- list(d0 = d0, d1 = d1, d0.ci = d0.ci, d1.ci = d1.ci, d0.p = d0.p, d1.p = d1.p, d0.sims = delta.0,d1.sims = delta.1,
                    z0 = z0, z1 = z1, z0.ci = z0.ci,z1.ci = z1.ci, z0.p = z0.p, z1.p = z1.p, z0.sims = zeta.0, z1.sims = zeta.1,
                    n0 = n0, n1 = n1, n0.ci = n0.ci, n1.ci = n1.ci, n0.p = n0.p, n1.p = n1.p, n0.sims = nu.0, n1.sims = nu.1,
                    tau.coef = tau.coef, tau.ci = tau.ci, tau.p = tau.p, tau.sims = tau,
                    d.avg = d.avg, d.avg.p = d.avg.p, d.avg.ci = d.avg.ci, d.avg.sims = delta.avg,
                    z.avg = z.avg, z.avg.p = z.avg.p, z.avg.ci = z.avg.ci, z.avg.sims = zeta.avg,
                    n.avg = n.avg, n.avg.p = n.avg.p, n.avg.ci = n.avg.ci, n.avg.sims = nu.avg,
                    treat = treat, mediator = mediator,
                    conf.level = conf.level,
                    model.y = model.y, model.m = lmodel.m,
                    control.value = control.value, treat.value = treat.value,
                    nobs = N, sims = J,cov=cov,cor=cor,
                    ORd0 = ORd0, ORd1 = ORd1, ORd0.ci = ORd0.ci, ORd1.ci = ORd1.ci, ORd0.p = ORd0.p, ORd1.p = ORd1.p, ORd0.sims = ORdelta.0,ORd1.sims = ORdelta.1,
                    ORz0 = ORz0, ORz1 = ORz1, ORz0.ci = ORz0.ci,ORz1.ci = ORz1.ci, ORz0.p = ORz0.p, ORz1.p = ORz1.p, ORz0.sims = ORzeta.0, ORz1.sims = ORzeta.1,
                    ORn0 = ORn0, ORn1 = ORn1, ORn0.ci = ORn0.ci, ORn1.ci = ORn1.ci, ORn0.p = ORn0.p, ORn1.p = ORn1.p, ORn0.sims = ORnu.0, ORn1.sims = ORnu.1,
                    ORtau.coef = ORtau.coef, ORtau.ci = ORtau.ci, ORtau.p = ORtau.p, ORtau.sims = ORtau,
                    ORd.avg = ORd.avg, ORd.avg.p = ORd.avg.p, ORd.avg.ci = ORd.avg.ci, ORd.avg.sims = ORdelta.avg,
                    ORz.avg = ORz.avg, ORz.avg.p = ORz.avg.p, ORz.avg.ci = ORz.avg.ci, ORz.avg.sims = ORzeta.avg,
                    ORn.avg = ORn.avg, ORn.avg.p = ORn.avg.p, ORn.avg.ci = ORn.avg.ci, ORn.avg.sims = ORnu.avg,
                    logORd0 = logORd0, logORd1 = logORd1, logORd0.ci = logORd0.ci, logORd1.ci = logORd1.ci, logORd0.p = logORd0.p, logORd1.p = logORd1.p, logORd0.sims = logORdelta.0,logORd1.sims = logORdelta.1,
                    logORz0 = logORz0, logORz1 = logORz1, logORz0.ci = logORz0.ci,logORz1.ci = logORz1.ci, logORz0.p = logORz0.p, logORz1.p = logORz1.p, logORz0.sims = logORzeta.0, logORz1.sims = logORzeta.1,
                    logORn0 = logORn0, logORn1 = logORn1, logORn0.ci = logORn0.ci, logORn1.ci = logORn1.ci, logORn0.p = logORn0.p, logORn1.p = logORn1.p, logORn0.sims = logORnu.0, logORn1.sims = logORnu.1,
                    logORtau.coef = logORtau.coef, logORtau.ci = logORtau.ci, logORtau.p = logORtau.p, logORtau.sims = logORtau,
                    logORd.avg = logORd.avg, logORd.avg.p = logORd.avg.p, logORd.avg.ci = logORd.avg.ci, logORd.avg.sims = logORdelta.avg,
                    logORz.avg = logORz.avg, logORz.avg.p = logORz.avg.p, logORz.avg.ci = logORz.avg.ci, logORz.avg.sims = logORzeta.avg,
                    logORn.avg = logORn.avg, logORn.avg.p = logORn.avg.p, logORn.avg.ci = logORn.avg.ci, logORn.avg.sims = logORnu.avg)
      }}
  }
  else {
    out <- list(d0 = d0, d1 = d1, d0.ci = d0.ci, d1.ci = d1.ci, d0.p = d0.p, d1.p = d1.p, d0.sims = delta.0,d1.sims = delta.1,
                d0.NM = d0.NM, d1.NM = d1.NM, d0.ci.NM = d0.ci.NM, d1.ci.NM = d1.ci.NM, d0.p.NM = d0.p.NM, d1.p.NM = d1.p.NM, d0.sims.NM = delta.0.NM, d1.sims.NM = delta.1.NM,
                z0 = z0, z1 = z1, z0.ci = z0.ci,z1.ci = z1.ci, z0.p = z0.p, z1.p = z1.p, z0.sims = zeta.0, z1.sims = zeta.1,
                n0 = n0, n1 = n1, n0.ci = n0.ci, n1.ci = n1.ci, n0.p = n0.p, n1.p = n1.p, n0.sims = nu.0, n1.sims = nu.1,
                n0.NM = n0.NM, n1.NM = n1.NM, n0.ci.NM = n0.ci.NM, n1.ci.NM = n1.ci.NM, n0.p.NM = n0.p.NM, n1.p.NM = n1.p.NM, n0.sims.NM = nu.0.NM, n1.sims.NM = nu.1.NM,
                tau.coef = tau.coef, tau.ci = tau.ci, tau.p = tau.p, tau.sims = tau,
                d.avg = d.avg, d.avg.p = d.avg.p, d.avg.ci = d.avg.ci, d.avg.sims = delta.avg,
                d.avg.NM = d.avg.NM, d.avg.p.NM = d.avg.p.NM, d.avg.ci.NM = d.avg.ci.NM, d.avg.sims.NM = delta.avg.NM,
                z.avg = z.avg, z.avg.p = z.avg.p, z.avg.ci = z.avg.ci, z.avg.sims = zeta.avg,
                n.avg = n.avg, n.avg.p = n.avg.p, n.avg.ci = n.avg.ci, n.avg.sims = nu.avg,
                n.avg.NM = n.avg.NM, n.avg.p.NM = n.avg.p.NM, n.avg.ci.NM = n.avg.ci.NM, n.avg.sims.NM = nu.avg.NM,
                treat = treat, mediator = mediator,
                conf.level = conf.level,
                model.y = model.y, model.m = lmodel.m,
                control.value = control.value, treat.value = treat.value,
                nobs = N, sims = J,cov=cov,cor=cor)

    if (!is.null(model.y$family)){
      if (model.y$family$link=="logit"){
        out <- list(d0 = d0, d1 = d1, d0.ci = d0.ci, d1.ci = d1.ci, d0.p = d0.p, d1.p = d1.p, d0.sims = delta.0,d1.sims = delta.1,
                    d0.NM = d0.NM, d1.NM = d1.NM, d0.ci.NM = d0.ci.NM, d1.ci.NM = d1.ci.NM, d0.p.NM = d0.p.NM, d1.p.NM = d1.p.NM, d0.sims.NM = delta.0.NM, d1.sims.NM = delta.1.NM,
                    z0 = z0, z1 = z1, z0.ci = z0.ci,z1.ci = z1.ci, z0.p = z0.p, z1.p = z1.p, z0.sims = zeta.0, z1.sims = zeta.1,
                    n0 = n0, n1 = n1, n0.ci = n0.ci, n1.ci = n1.ci, n0.p = n0.p, n1.p = n1.p, n0.sims = nu.0, n1.sims = nu.1,
                    n0.NM = n0.NM, n1.NM = n1.NM, n0.ci.NM = n0.ci.NM, n1.ci.NM = n1.ci.NM, n0.p.NM = n0.p.NM, n1.p.NM = n1.p.NM, n0.sims.NM = nu.0.NM, n1.sims.NM = nu.1.NM,
                    tau.coef = tau.coef, tau.ci = tau.ci, tau.p = tau.p, tau.sims = tau,
                    d.avg = d.avg, d.avg.p = d.avg.p, d.avg.ci = d.avg.ci, d.avg.sims = delta.avg,
                    d.avg.NM = d.avg.NM, d.avg.p.NM = d.avg.p.NM, d.avg.ci.NM = d.avg.ci.NM, d.avg.sims.NM = delta.avg.NM,
                    z.avg = z.avg, z.avg.p = z.avg.p, z.avg.ci = z.avg.ci, z.avg.sims = zeta.avg,
                    n.avg = n.avg, n.avg.p = n.avg.p, n.avg.ci = n.avg.ci, n.avg.sims = nu.avg,
                    n.avg.NM = n.avg.NM, n.avg.p.NM = n.avg.p.NM, n.avg.ci.NM = n.avg.ci.NM, n.avg.sims.NM = nu.avg.NM,
                    treat = treat, mediator = mediator,
                    conf.level = conf.level,
                    model.y = model.y, model.m = lmodel.m,
                    control.value = control.value, treat.value = treat.value,
                    nobs = N, sims = J,cov=cov,cor=cor,
                    ORd0 = ORd0, ORd1 = ORd1, ORd0.ci = ORd0.ci, ORd1.ci = ORd1.ci, ORd0.p = ORd0.p, ORd1.p = ORd1.p, ORd0.sims = ORdelta.0,ORd1.sims = ORdelta.1,
                    ORd0.NM = ORd0.NM, ORd1.NM = ORd1.NM, ORd0.ci.NM = ORd0.ci.NM, ORd1.ci.NM = ORd1.ci.NM, ORd0.p.NM = ORd0.p.NM, ORd1.p.NM = ORd1.p.NM, ORd0.sims.NM = ORdelta.0.NM, ORd1.sims.NM = ORdelta.1.NM,
                    ORz0 = ORz0, ORz1 = ORz1, ORz0.ci = ORz0.ci,ORz1.ci = ORz1.ci, ORz0.p = ORz0.p, ORz1.p = ORz1.p, ORz0.sims = ORzeta.0, ORz1.sims = ORzeta.1,
                    ORn0 = ORn0, ORn1 = ORn1, ORn0.ci = ORn0.ci, ORn1.ci = ORn1.ci, ORn0.p = ORn0.p, ORn1.p = ORn1.p, ORn0.sims = ORnu.0, ORn1.sims = ORnu.1,
                    ORn0.NM = ORn0.NM, ORn1.NM = ORn1.NM, ORn0.ci.NM = ORn0.ci.NM, ORn1.ci.NM = ORn1.ci.NM, ORn0.p.NM = ORn0.p.NM, ORn1.p.NM = ORn1.p.NM, ORn0.sims.NM = ORnu.0.NM, ORn1.sims.NM = ORnu.1.NM,
                    ORtau.coef = ORtau.coef, ORtau.ci = ORtau.ci, ORtau.p = ORtau.p, ORtau.sims = ORtau,
                    ORd.avg = ORd.avg, ORd.avg.p = ORd.avg.p, ORd.avg.ci = ORd.avg.ci, ORd.avg.sims = ORdelta.avg,
                    ORd.avg.NM = ORd.avg.NM, ORd.avg.p.NM = ORd.avg.p.NM, ORd.avg.ci.NM = ORd.avg.ci.NM, ORd.avg.sims.NM = ORdelta.avg.NM,
                    ORz.avg = ORz.avg, ORz.avg.p = ORz.avg.p, ORz.avg.ci = ORz.avg.ci, ORz.avg.sims = ORzeta.avg,
                    ORn.avg = ORn.avg, ORn.avg.p = ORn.avg.p, ORn.avg.ci = ORn.avg.ci, ORn.avg.sims = ORnu.avg,
                    ORn.avg.NM = ORn.avg.NM, ORn.avg.p.NM = ORn.avg.p.NM, ORn.avg.ci.NM = ORn.avg.ci.NM, ORn.avg.sims.NM = ORnu.avg.NM,
                    logORd0 = logORd0, logORd1 = logORd1, logORd0.ci = logORd0.ci, logORd1.ci = logORd1.ci, logORd0.p = logORd0.p, logORd1.p = logORd1.p, logORd0.sims = logORdelta.0,logORd1.sims = logORdelta.1,
                    logORd0.NM = logORd0.NM, logORd1.NM = logORd1.NM, logORd0.ci.NM = logORd0.ci.NM, logORd1.ci.NM = logORd1.ci.NM, logORd0.p.NM = logORd0.p.NM, logORd1.p.NM = logORd1.p.NM, logORd0.sims.NM = logORdelta.0.NM, logORd1.sims.NM = logORdelta.1.NM,
                    logORz0 = logORz0, logORz1 = logORz1, logORz0.ci = logORz0.ci,logORz1.ci = logORz1.ci, logORz0.p = logORz0.p, logORz1.p = logORz1.p, logORz0.sims = logORzeta.0, logORz1.sims = logORzeta.1,
                    logORn0 = logORn0, logORn1 = logORn1, logORn0.ci = logORn0.ci, logORn1.ci = logORn1.ci, logORn0.p = logORn0.p, logORn1.p = logORn1.p, logORn0.sims = logORnu.0, logORn1.sims = logORnu.1,
                    logORn0.NM = logORn0.NM, logORn1.NM = logORn1.NM, logORn0.ci.NM = logORn0.ci.NM, logORn1.ci.NM = logORn1.ci.NM, logORn0.p.NM = logORn0.p.NM, logORn1.p.NM = logORn1.p.NM, logORn0.sims.NM = logORnu.0.NM, logORn1.sims.NM = logORnu.1.NM,
                    logORtau.coef = logORtau.coef, logORtau.ci = logORtau.ci, logORtau.p = logORtau.p, logORtau.sims = logORtau,
                    logORd.avg = logORd.avg, logORd.avg.p = logORd.avg.p, logORd.avg.ci = logORd.avg.ci, logORd.avg.sims = logORdelta.avg,
                    logORd.avg.NM = logORd.avg.NM, logORd.avg.p.NM = logORd.avg.p.NM, logORd.avg.ci.NM = logORd.avg.ci.NM, logORd.avg.sims.NM = logORdelta.avg.NM,
                    logORz.avg = logORz.avg, logORz.avg.p = logORz.avg.p, logORz.avg.ci = logORz.avg.ci, logORz.avg.sims = logORzeta.avg,
                    logORn.avg = logORn.avg, logORn.avg.p = logORn.avg.p, logORn.avg.ci = logORn.avg.ci, logORn.avg.sims = logORnu.avg,
                    logORn.avg.NM = logORn.avg.NM, logORn.avg.p.NM = logORn.avg.p.NM, logORn.avg.ci.NM = logORn.avg.ci.NM, logORn.avg.sims.NM = logORnu.avg.NM)
      }}
  }
  class(out) <- "mm"
  return(out)
}
