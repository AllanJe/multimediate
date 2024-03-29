#' blockmediate
#'
#'
#' 'blockmediate' is used to estimate various quantities for causal mediation analysis as "multimediate", not for several mediators but for for different block of mediators.
#'
#'@param lmodel.m list of fitted models object for mediators. Can be of class 'lm', 'polr','glm'.
#'@param model.y a fitted model object for the outcome. Can be of class 'lm', 'polr','glm'.
#'@param treat a character string indicating the name of the treatment variable used in the models. The treatment can be either binary (integer or a two-valued factor) or continuous (numeric).
#'@param treat.value value of the treatment variable used as the treatment condition. Default is 1.
#'@param control.value value of the treatment variable used as the control condition. Default is 0.
#'@param J number of Monte Carlo draws for quasi-Bayesian approximation.
#'@param conf.level level of the returned two-sided confidence intervals. Default is to return the 2.5 and 97.5 percentiles of the simulated quantities.
#'@param clust vector of length equal to the number of mediator. The order of the numbers corresponds to the order of the mediators in lmodel.m and the mediators belonging to the same cluster have the same number. Default is 'NULL', clusters are define using clustering spectral.
#'@param seuil threshold correlation value for spectral clustering. Default is 0.
#'
#'
#'@return blockmediate returns an object of class "bm", a list that contains the components listed below.
#' The function summary (i.e., summary.bm) can be used to obtain a table of the results.
#' @export
#'
#' @importFrom MASS mvrnorm
#' @importFrom mvtnorm rmvnorm

blockmediate = function(lmodel.m,model.y,treat,treat.value=1,control.value=0,J=1000,conf.level=0.95,clust=NULL,seuil=0){
  N=dim(lmodel.m[[1]]$model)[1]
  NUM=length(lmodel.m)

  mediator=c()
  MModel = list()
  for (num in 1:NUM){
    MModel[[num]]=rmvnorm(J, mean = c(coef(lmodel.m[[num]]),lmodel.m[[num]]$zeta), sigma = vcov(lmodel.m[[num]]))
    mediator=c(mediator,names(lmodel.m[[num]]$model)[1])
  }

  mcov=CorCond(e=10^(-10),lmodel.m)

  Corre=mcov$correstim
  if(is.null(clust)){
    #clust=SpecClust(mcov$correstim,seuil)
    for(s in 1:NUM){
      for(d in 1:NUM){
        if(abs(Corre[s,d])>=0.99){
          Corre[s,d]=(-1/2)*log(1-0.99^2)}
        else{
          Corre[s,d]=(-1/2)*log(1-Corre[s,d]^2)
        }
      }}

    D <- diag(apply(Corre, 1, sum))
    UL <- D - Corre
    round(UL,1)
    evL <- eigen(UL, symmetric=TRUE)
    par(mfrow=c(1,1))
    barplot(evL$values)
    k=as.numeric(readline(prompt="How many cluster?"))
    Z   <- evL$vectors[,(ncol(evL$vectors)-k+1):ncol(evL$vectors)]
    clust=kmeans(Z,k,nstart=3)$cluster
  }
  NM=max(clust)
  Sigma=Sigma_adjusted=mcov$sigmaestim

  for(i in 1:NUM){
    for(j in 1:NUM){
      if (clust[i]!=clust[j]){
        Sigma_adjusted[i,j]=0
      }
    }
  }
  error=mvrnorm(n=N*J,mu=rep(0,NUM),Sigma=Sigma_adjusted,tol=100)



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
  PredictM1<-PredictM0<-PredictM1b<-PredictM0b<- array(0, dim=c(J,N,NUM))


  for (num in 1:NUM){
    pred.data.t <- pred.data.c <- model.frame(lmodel.m[[num]])

    if (is.factor(lmodel.m[[num]]$model[,treat])) {
      pred.data.t[, treat] <- factor(treat.value, levels = levels(lmodel.m[[num]]$model[, treat]))
      pred.data.c[, treat] <- factor(control.value, levels = levels(lmodel.m[[num]]$model[, treat]))
    }
    else {
      pred.data.t[, treat] <- treat.value
      pred.data.c[, treat] <- control.value
    }


    if(inherits(lmodel.m[[num]],"polr")){
      mmat.t <- model.matrix(terms(lmodel.m[[num]]), data = pred.data.t)
      mmat.c <- model.matrix(terms(lmodel.m[[num]]), data = pred.data.c)
      mmat.t=mmat.t[,-1]
      mmat.c=mmat.c[,-1]

      if (is.null(dim(mmat.t))){
        muM1 <- tcrossprod(MModel[[num]][,1], mmat.t)
        muM0 <- tcrossprod(MModel[[num]][,1], mmat.c)
      }
      else{
        muM1 <- tcrossprod(MModel[[num]][,1:dim(mmat.t)[2]], mmat.t)
        muM0 <- tcrossprod(MModel[[num]][,1:dim(mmat.t)[2]], mmat.c)
      }


      PredictM1[,,num] <- array(muM1,dim=c(J,N)) + array(error[,num], dim=c(J,N))
      PredictM0[,,num] <- array(muM0,dim=c(J,N)) + array(error[,num], dim=c(J,N))
      PredictM1b[,,num] = array(muM1,dim=c(J,N)) + array(error[,num], dim=c(J,N))
      PredictM0b[,,num] = array(muM0,dim=c(J,N)) + array(error[,num], dim=c(J,N))
      if (is.null(dim(mmat.t))){
        seuil=cbind(-Inf,MModel[[num]][,-1],Inf)}
      else{
        seuil=cbind(-Inf,MModel[[num]][,-(1:dim(mmat.t)[2])],Inf)
      }
      for (k in 1:length(lmodel.m[[num]]$lev)){
        for (n in 1:N){
          a=which(PredictM1b[,n,num]>seuil[,k] & PredictM1b[,n,num]<=seuil[,k+1])
          b=which(PredictM0b[,n,num]>seuil[,k] & PredictM0b[,n,num]<=seuil[,k+1])
          PredictM1[a,n,num]=lmodel.m[[num]]$lev[k]
          PredictM0[b,n,num]=lmodel.m[[num]]$lev[k]
        }
      }
    }
    else{
      mmat.t <- model.matrix(terms(lmodel.m[[num]]), data = pred.data.t)
      mmat.c <- model.matrix(terms(lmodel.m[[num]]), data = pred.data.c)

      muM1 <- tcrossprod(MModel[[num]], mmat.t)
      muM0 <- tcrossprod(MModel[[num]], mmat.c)

      PredictM1[,,num] <- array(muM1,dim=c(J,N)) + array(error[,num], dim=c(J,N))
      PredictM0[,,num] <- array(muM0,dim=c(J,N)) + array(error[,num], dim=c(J,N))
      if (inherits(lmodel.m[[num]],"glm")){
        PredictM1[,,num]=(PredictM1[,,num]>0)*1
        PredictM0[,,num]=(PredictM0[,,num]>0)*1
      }
    }


  }

  effect.tmp.NM=array(NA, dim = c(N, J, 2, NM))
  effect.tmp=array(NA, dim = c(N, J, 4))
  OR.NM=array(NA, dim = c(J, 2, NM))
  OR=array(NA, dim = c(J, 4))

  for (e in 1:4) {
    tt <- switch(e, c(1, 1, 1, 0), c(0, 0, 1, 0),
                 c(1, 0, 1, 1), c(1, 0, 0, 0))
    Pr0<-Pr1<-ORPr0<-ORPr1<- matrix(NA, nrow = N, ncol = J)
    if (NM!=1 & e<=2) {
      Pr0.NM<-Pr1.NM <-ORPr0.NM<-ORPr1.NM<- array(NA,dim=c(N,J,NM))
    }

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




      for (num in 1:NUM){
        if(inherits(lmodel.m[[num]],"polr")){
          pred.data.c[,mediator[num]]=as.factor(pred.data.c[,mediator[num]])
          levels(pred.data.c[,mediator[num]])=lmodel.m[[num]]$lev
          pred.data.t[,mediator[num]]=as.factor(pred.data.t[,mediator[num]])
          levels(pred.data.t[,mediator[num]])=lmodel.m[[num]]$lev
        }

        if(inherits(lmodel.m[[num]],"lm") & !inherits(lmodel.m[[num]],"glm")){
          pred.data.c[,mediator[num]]=as.numeric(pred.data.c[,mediator[num]])
          pred.data.t[,mediator[num]]=as.numeric(pred.data.t[,mediator[num]])
        }

        if(inherits(lmodel.m[[num]],"glm")){
          pred.data.c[,mediator[num]]=as.factor(pred.data.c[,mediator[num]])
          levels(pred.data.c[,mediator[num]])=levels(lmodel.m[[num]]$model[,1])
          pred.data.t[,mediator[num]]=as.factor(pred.data.t[,mediator[num]])
          levels(pred.data.t[,mediator[num]])=levels(lmodel.m[[num]]$model[,1])
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
      }
      else{
        ymat.t=ymat.t[,-1]
        ymat.c=ymat.c[,-1]
        Pr1[,j] <- t(as.matrix(YModel[j,1:dim(ymat.t)[2]])) %*% t(ymat.t)+ Yerror
        Pr0[,j] <- t(as.matrix(YModel[j,1:dim(ymat.c)[2]])) %*% t(ymat.c)+ Yerror

        seuilY=cbind(-Inf,YModel[,-(1:dim(ymat.t)[2]+1)],Inf)
        Pr1b=as.numeric(Pr1[,j])
        Pr0b=as.numeric(Pr0[,j])
        for (k in 1:length(model.y$lev)){
          a=which(Pr1b>seuilY[j,k] & Pr1b<=seuilY[j,k+1])
          b=which(Pr0b>seuilY[j,k] & Pr0b<=seuilY[j,k+1])
          Pr1[a,j]=model.y$lev[k]
          Pr0[b,j]=model.y$lev[k]
        }
        Pr1=matrix(as.numeric(Pr1),N,J)
        Pr0=matrix(as.numeric(Pr0),N,J)
      }

      if (NM!=1 & e<=2) {

        for (nm in 1:NM){

          inm=which(clust==nm)
          if (tt[3]==1){
            PredictMt=PredictM1[j,,inm]
          }
          else {
            PredictMt=PredictM0[j,,inm]
          }
          if (tt[1]==1){
            PredictWt=PredictM1[j,,-inm]
          }
          else {
            PredictWt=PredictM0[j,,-inm]
          }

          if (tt[4]==1){
            PredictMc=PredictM1[j,,inm]
          }
          else {
            PredictMc=PredictM0[j,,inm]
          }
          if (tt[2]==1){
            PredictWc=PredictM1[j,,-inm]
          }
          else {
            PredictWc=PredictM0[j,,-inm]
          }


          pred.data.t[, mediator[inm]] <- PredictMt
          pred.data.t[, mediator[-inm]] <- PredictWt
          pred.data.c[, mediator[inm]] <- PredictMc
          pred.data.c[, mediator[-inm]] <- PredictWc
          for (numt in 1:NUM){
            if(inherits(lmodel.m[[numt]],"polr")){

              pred.data.c[,mediator[numt]]=as.factor(pred.data.c[,mediator[numt]])
              levels(pred.data.c[,mediator[numt]])=lmodel.m[[numt]]$lev
              pred.data.t[,mediator[numt]]=as.factor(pred.data.t[,mediator[numt]])
              levels(pred.data.t[,mediator[numt]])=lmodel.m[[numt]]$lev
            }

            if(inherits(lmodel.m[[numt]],"lm") & !inherits(lmodel.m[[numt]],"glm")){
              pred.data.c[,mediator[numt]]=as.numeric(pred.data.c[,mediator[numt]])
              pred.data.t[,mediator[numt]]=as.numeric(pred.data.t[,mediator[numt]])
            }

            if(inherits(lmodel.m[[numt]],"glm")){
              pred.data.c[,mediator[numt]]=as.factor(pred.data.c[,mediator[numt]])
              levels(pred.data.c[,mediator[numt]])=levels(lmodel.m[[numt]]$model[,1])
              pred.data.t[,mediator[numt]]=as.factor(pred.data.t[,mediator[numt]])
              levels(pred.data.t[,mediator[numt]])=levels(lmodel.m[[numt]]$model[,1])
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
          }
          else{
            ymat.t=ymat.t[,-1]
            ymat.c=ymat.c[,-1]
            Pr1.NM[,j,nm] <- t(as.matrix(YModel[j,1:dim(ymat.t)[2]])) %*% t(ymat.t)+ Yerror
            Pr0.NM[,j,nm] <- t(as.matrix(YModel[j,1:dim(ymat.c)[2]])) %*% t(ymat.c)+ Yerror

            Pr1b.NM=as.numeric(Pr1.NM[,j,nm])
            Pr0b.NM=as.numeric(Pr0.NM[,j,nm])
            for (k in 1:length(model.y$lev)){
              a=which(Pr1b.NM>seuilY[j,k] & Pr1b.NM<=seuilY[j,k+1])
              b=which(Pr0b.NM>seuilY[j,k] & Pr0b.NM<=seuilY[j,k+1])
              Pr1.NM[a,j,nm]=model.y$lev[k]
              Pr0.NM[b,j,nm]=model.y$lev[k]
            }
            Pr1.NM=array(as.numeric(Pr1.NM),dim=c(N,J,NM))
            Pr0.NM=array(as.numeric(Pr0.NM),dim=c(N,J,NM))
          }
        }
      }


    }

    if (!is.null(model.y$family)){

      effect.tmp[,,e]=(Pr1>0)-(Pr0>0)
      if (model.y$family$link=="logit"){
        OR[,e]=(1-apply(Pr0>0,c(2),mean))/(1-apply(Pr1>0,c(2),mean))*(1+(apply(effect.tmp[,,e],c(2),mean)/apply(Pr0>0,c(2),mean)))
      }
      if (NM!=1 & e<=2){

        effect.tmp.NM[,,e,]=(Pr1.NM>0)-(Pr0.NM>0)
        if (model.y$family$link=="logit"){
          OR.NM[,e,]=(1-apply(Pr0.NM>=0, c(2,3), mean))/(1-apply(Pr1.NM>=0, c(2,3), mean))*(1+(apply(effect.tmp.NM[,,e,],c(2,3),mean)/apply(Pr0.NM>=0, c(2,3), mean)))
        }
      }
    }

    else {
      effect.tmp[,,e]=Pr1-Pr0
      if (NM!=1 & e<=2){
        effect.tmp.NM[,,e,]=Pr1.NM-Pr0.NM
      }
    }
  }



  et1 <- effect.tmp[, , 1]
  et2 <- effect.tmp[, , 2]
  et3 <- effect.tmp[, , 3]
  et4 <- effect.tmp[, , 4]
  et5 <- effect.tmp.NM[, ,1,]
  et6 <- effect.tmp.NM[, ,2,]

  if (!is.null(model.y$family)){
    if ( model.y$family$link=="logit"){


      ORdelta.1 <- OR[ , 1]
      ORdelta.0 <- OR[ , 2]
      ORzeta.1 <- OR[ , 3]
      ORzeta.0 <- OR[ , 4]
      ORdelta.1.NM <- OR.NM[ ,1,]
      ORdelta.0.NM <- OR.NM[ ,2,]

    }}


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
  d0 <- mean(delta.0,na.rm=TRUE)
  d1 <- mean(delta.1,na.rm=TRUE)
  z1 <- mean(zeta.1,na.rm=TRUE)
  z0 <- mean(zeta.0,na.rm=TRUE)
  tau.coef <- mean(tau,na.rm=TRUE)
  n0 <- median(nu.0,na.rm=TRUE)
  n1 <- median(nu.1,na.rm=TRUE)
  d.avg <- (d0 + d1)/2
  z.avg <- (z0 + z1)/2
  n.avg <- (n0 + n1)/2

  if (!is.null(model.y$family)){
    if ( model.y$family$link=="logit"){
      ORtau <- (ORzeta.1*ORdelta.0 + ORzeta.0*ORdelta.1)/2
      ORnu.1 <- ORdelta.1/ORtau
      ORnu.0 <- ORdelta.0/ORtau


      ORtau.coef <- mean(ORtau,na.rm=TRUE)


      ORdelta.avg <- (ORdelta.1 + ORdelta.0)/2
      ORnu.avg <- (ORnu.1 + ORnu.0)/2
      ORzeta.avg <- (ORzeta.1 + ORzeta.0)/2

      ORd1 <- mean(ORdelta.1,na.rm=TRUE)
      ORn1 <- median(ORnu.1,na.rm=TRUE)
      ORd0 <- mean(ORdelta.0,na.rm=TRUE)
      ORn0 <- median(ORnu.0,na.rm=TRUE)
      ORz1 <- mean(ORzeta.1,na.rm=TRUE)
      ORz0 <- mean(ORzeta.0,na.rm=TRUE)

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

      logORtau.coef <- mean(logORtau,na.rm=TRUE)
      logORnu.avg <- (logORnu.1 + logORnu.0)/2
      logORn0 <- median(logORnu.0,na.rm=TRUE)
      logORn1 <- median(logORnu.1,na.rm=TRUE)
      logORn.avg <- (logORn0 + logORn1)/2

      logORdelta.avg <- (logORdelta.1 + logORdelta.0)/2
      logORzeta.avg <- (logORzeta.1 + logORzeta.0)/2
      logORd0 <- mean(logORdelta.0,na.rm=TRUE)
      logORd1 <- mean(logORdelta.1,na.rm=TRUE)
      logORz1 <- mean(logORzeta.1,na.rm=TRUE)
      logORz0 <- mean(logORzeta.0,na.rm=TRUE)
      logORtau.coef <- mean(logORtau,na.rm=TRUE)
      logORd.avg <- (logORd0 + logORd1)/2
      logORz.avg <- (logORz0 + logORz1)/2
    }}

  if (NM!=1){
    delta.1.NM <- t((apply(et5, c(2,3), mean)))
    delta.0.NM <- t(as.matrix(apply(et6, c(2,3), mean)))

    nu.0.NM <- delta.0.NM/array(tau,dim=c(NM,J))
    nu.1.NM <- delta.1.NM/array(tau,dim=c(NM,J))
    delta.avg.NM <- (delta.1.NM + delta.0.NM)/2
    nu.avg.NM <- (nu.1.NM + nu.0.NM)/2
    d0.NM <- apply(delta.0.NM,1,mean,na.rm=TRUE)
    d1.NM <- apply(delta.1.NM,1,mean,na.rm=TRUE)
    n0.NM <- apply(nu.0.NM,1,median,na.rm=TRUE)
    n1.NM <- apply(nu.1.NM,1,median,na.rm=TRUE)
    d.avg.NM <- (d0.NM + d1.NM)/2
    n.avg.NM <- (n0.NM + n1.NM)/2

    if (!is.null(model.y$family)){
      if (model.y$family$link=="logit"){

        ORnu.0.NM <- ORdelta.0.NM/array(ORtau,dim=c(J,NM))
        ORnu.1.NM <- ORdelta.1.NM/array(ORtau,dim=c(J,NM))
        ORdelta.avg.NM <- (ORdelta.1.NM + ORdelta.0.NM)/2
        ORnu.avg.NM <- (ORnu.1.NM + ORnu.0.NM)/2
        ORd0.NM <- apply(ORdelta.0.NM,2,mean,na.rm=TRUE)
        ORd1.NM <- apply(ORdelta.1.NM,2,mean,na.rm=TRUE)
        ORn0.NM <- apply(ORnu.0.NM,2,median,na.rm=TRUE)
        ORn1.NM <- apply(ORnu.1.NM,2,median,na.rm=TRUE)
        ORd.avg.NM <- (ORd0.NM + ORd1.NM)/2
        ORn.avg.NM <- (ORn0.NM + ORn1.NM)/2

        logORdelta.1.NM <- log(ORdelta.1.NM)
        logORdelta.0.NM <- log(ORdelta.0.NM)
        logORnu.0.NM <- logORdelta.0.NM/array(logORtau,dim=c(J,NM))
        logORnu.1.NM <- logORdelta.1.NM/array(logORtau,dim=c(J,NM))

        logORnu.avg.NM <- (logORnu.1.NM + logORnu.0.NM)/2
        logORn0.NM <- apply(logORdelta.0.NM/logORtau.coef,2,median,na.rm=TRUE)
        logORn1.NM <- apply(logORdelta.1.NM/logORtau.coef,2,median,na.rm=TRUE)
        logORn.avg.NM <- (logORn0.NM + logORn1.NM)/2

        logORdelta.avg.NM <- (logORdelta.1.NM + logORdelta.0.NM)/2
        logORd0.NM <- apply(logORdelta.0.NM,2,mean,na.rm=TRUE)
        logORd1.NM <- apply(logORdelta.1.NM,2,mean,na.rm=TRUE)
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

      ORd0.p <- pval(ORdelta.0, ORd0)
      ORd1.p <- pval(ORdelta.1, ORd1)
      ORd.avg.p <- pval(ORdelta.avg, ORd.avg)
      ORz0.p <- pval(ORzeta.0, ORz0)
      ORz1.p <- pval(ORzeta.1, ORz1)
      ORz.avg.p <- pval(ORzeta.avg, ORz.avg)
      ORn0.p <- pval(ORnu.0, ORn0)
      ORn1.p <- pval(ORnu.1, ORn1)
      ORn.avg.p <- pval(ORnu.avg, ORn.avg)
      ORtau.p <- pval(ORtau, ORtau.coef)

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
    for (num in 1:NM){
      d0.ci.NM[num,] <- quantile(delta.0.NM[num,], c(low, high), na.rm = TRUE)
      d1.ci.NM[num,] <- quantile(delta.1.NM[num,], c(low, high), na.rm = TRUE)
      n0.ci.NM[num,] <- quantile(nu.0.NM[num,], c(low, high), na.rm = TRUE)
      n1.ci.NM[num,] <- quantile(nu.1.NM[num,], c(low, high), na.rm = TRUE)
      d.avg.ci.NM[num,] <- quantile(delta.avg.NM[num,], c(low, high), na.rm = TRUE)
      n.avg.ci.NM[num,] <- quantile(nu.avg.NM[num,], c(low, high), na.rm = TRUE)

      d0.p.NM[num] <- pval(delta.0.NM[num,], d0.NM[num])
      d1.p.NM[num] <- pval(delta.1.NM[num,], d1.NM[num])
      d.avg.p.NM[num] <- pval(delta.avg.NM[num,], d.avg.NM[num])
      n0.p.NM[num] <- pval(nu.0.NM[num,], n0.NM[num])
      n1.p.NM[num] <- pval(nu.1.NM[num,], n1.NM[num])
      n.avg.p.NM[num] <- pval(nu.avg.NM[num,], n.avg.NM[num])
    }

    if (!is.null(model.y$family)){
      if (model.y$family$link=="logit"){
        ORd0.ci.NM <- ORd1.ci.NM <- ORn0.ci.NM <- ORn1.ci.NM <-ORd.avg.ci.NM <- ORn.avg.ci.NM <- array(NA,dim=c(NM,2))
        ORd0.p.NM <- ORd1.p.NM <- ORd.avg.p.NM <- ORn0.p.NM <- ORn1.p.NM <- ORn.avg.p.NM<- rep(NA,NM)

        logORd0.ci.NM <- logORd1.ci.NM <- logORn0.ci.NM <- logORn1.ci.NM <- logORd.avg.ci.NM <- logORn.avg.ci.NM <- array(NA,dim=c(NM,2))
        logORd0.p.NM <- logORd1.p.NM <- logORd.avg.p.NM <- logORn0.p.NM <- logORn1.p.NM <- logORn.avg.p.NM <- rep(NA,NM)
        for (num in 1:NM){
          ORd0.ci.NM[num,] <- quantile(ORdelta.0.NM[num,], c(low, high), na.rm = TRUE)
          ORd1.ci.NM[num,] <- quantile(ORdelta.1.NM[num,], c(low, high), na.rm = TRUE)
          ORn0.ci.NM[num,] <- quantile(ORnu.0.NM[num,], c(low, high), na.rm = TRUE)
          ORn1.ci.NM[num,] <- quantile(ORnu.1.NM[num,], c(low, high), na.rm = TRUE)
          ORd.avg.ci.NM[num,] <- quantile(ORdelta.avg.NM[num,], c(low, high), na.rm = TRUE)
          ORn.avg.ci.NM[num,] <- quantile(ORnu.avg.NM[num,], c(low, high), na.rm = TRUE)

          ORd0.p.NM[num] <- pval(ORdelta.0.NM[num,], ORd0.NM[num])
          ORd1.p.NM[num] <- pval(ORdelta.1.NM[num,], ORd1.NM[num])
          ORd.avg.p.NM[num] <- pval(ORdelta.avg.NM[num,], ORd.avg.NM[num])
          ORn0.p.NM[num] <- pval(ORnu.0.NM[num,], ORn0.NM[num])
          ORn1.p.NM[num] <- pval(ORnu.1.NM[num,], ORn1.NM[num])
          ORn.avg.p.NM[num] <- pval(ORnu.avg.NM[num,], ORn.avg.NM[num])

          logORd0.ci.NM[num,] <- quantile(logORdelta.0.NM[num,], c(low, high), na.rm = TRUE)
          logORd1.ci.NM[num,] <- quantile(logORdelta.1.NM[num,], c(low, high), na.rm = TRUE)
          logORd.avg.ci.NM[num,] <- quantile(logORdelta.avg.NM[num,], c(low, high), na.rm = TRUE)
          logORn0.ci.NM[num,] <- quantile(logORnu.0.NM[num,], c(low, high), na.rm = TRUE)
          logORn1.ci.NM[num,] <- quantile(logORnu.1.NM[num,], c(low, high), na.rm = TRUE)
          logORn.avg.ci.NM[num,] <- quantile(logORnu.avg.NM[num,], c(low, high), na.rm = TRUE)
          logORd0.p.NM[num] <- pval(logORdelta.0.NM[num,], logORd0.NM[num])
          logORd1.p.NM[num] <- pval(ORdelta.1.NM[num,], logORd1.NM[num])
          logORd.avg.p.NM[num] <- pval(logORdelta.avg.NM[num,], logORd.avg.NM[num])
          logORn0.p.NM[num] <- pval(logORnu.0.NM[num,], logORn0.NM[num])
          logORn1.p.NM[num] <- pval(logORnu.1.NM[num,], logORn1.NM[num])
          logORn.avg.p.NM[num] <- pval(logORnu.avg.NM[num,], logORn.avg.NM[num])
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
                nobs = N, sims = J,clust=clust,sigma=Sigma,sigma_adjusted=Sigma_adjusted,corr=mcov$correstim)
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
                    nobs = N, sims = J,
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
                    logORn.avg = logORn.avg, logORn.avg.p = logORn.avg.p, logORn.avg.ci = logORn.avg.ci, logORn.avg.sims = logORnu.avg,clust=clust,sigma=Sigma,sigma_adjusted=Sigma_adjusted,corr=mcov$correstim)
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
                nobs = N, sims = J,clust=clust,sigma=Sigma,sigma_adjusted=Sigma_adjusted,corr=mcov$correstim)

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
                    nobs = N, sims = J,
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
                    logORn.avg.NM = logORn.avg.NM, logORn.avg.p.NM = logORn.avg.p.NM, logORn.avg.ci.NM = logORn.avg.ci.NM, logORn.avg.sims.NM = logORnu.avg.NM,clust=clust,sigma=Sigma,sigma_adjusted=Sigma_adjusted,corr=mcov$correstim)
      }}
  }
  class(out) <- "bm"
  return(out)
}
