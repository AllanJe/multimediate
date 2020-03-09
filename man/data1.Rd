\name{data1}
\alias{data1}
\docType{data}
\title{
Toy dataset 1
}

\description{A toy dataset for causal mediation analysis with 3 uncausally related mediators.
}

\usage{data("data1")}
\format{
  A data frame with 1000 observations on the following 8 variables.
  \describe{
    \item{\code{Treatment}}{a factor with levels \code{0} \code{1}}
    \item{\code{C1}}{a factor with levels \code{0} \code{1}}
    \item{\code{C2}}{a factor with levels \code{0} \code{1}}
    \item{\code{C3}}{a factor with levels \code{0} \code{1}}
    \item{\code{M1}}{a numeric vector}
    \item{\code{M2}}{a numeric vector}
    \item{\code{M3}}{a numeric vector}
    \item{\code{Outcome}}{a numeric vector}
  }
}

\examples{
#data(data1)
#data1$Treatment=as.factor(data1$Treatment)
#data1$C1=as.factor(data1$C1)
#data1$C2=as.factor(data1$C2)
#data1$C3=as.factor(data1$C3)

#M1reg=lm(M1~ Treatment + C1 + C2, data=data1)
#M2reg=lm(M2~ Treatment + C1 + C3, data=data1)
#M3reg=lm(M3~ Treatment + C1     , data=data1)

#Yreg=lm(Outcome~ Treatment + M1 + M2 + M3 + C1 + C2 + C3, data=data1)

#med.analysis=multimediate(lmodel.m=list(M1reg,M2reg,M3reg),correlated=TRUE,
                          #model.y=Yreg,treat="Treatment",treat.value=1,
                          #control.value=0,J=1000,conf.level=0.95)

#summary.mm(out=med.analysis,opt="avg")
}
\keyword{datasets}
