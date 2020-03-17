\name{data2}
\alias{data2}
\docType{data}
\title{
Toy dataset 2
}

\description{A toy dataset for causal mediation analysis with 3 uncausally related mediators.
}
\usage{data("data2")}
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
    \item{\code{Outcome}}{a factor with levels \code{FALSE} \code{TRUE}}
  }
}

\examples{
#data(data2)
#data2$Treatment=as.factor(data2$Treatment)
#data2$C1=as.factor(data2$C1)
#data2$C2=as.factor(data2$C2)
#data2$C3=as.factor(data2$C3)
#data2$M1=as.numeric(data2$M1)
#data2$M2=as.numeric(data2$M2)
#data2$M3=as.numeric(data2$M3)
#data2$Outcome=as.factor(data2$Outcome)
#summary(data2)

#M1reg=lm(M1~ Treatment + C1, data=data2)
#M2reg=lm(M2~ Treatment + C2, data=data2)
#M3reg=lm(M3~ Treatment + C3, data=data2)

#Yreg=glm(Outcome~ Treatment + M1 + M2 + M3 + C1 + C2 + C3, data=data2, family = binomial("logit"))

#med.analysis=multimediate(lmodel.m=list(M1reg,M2reg,M3reg),
#correlated=TRUE,model.y=Yreg,treat="Treatment",treat.value=1,control.value=0,J=1000,conf.level=0.95)

# summary(med.analysis,opt="avg",logit="effects")
# summary(med.analysis,opt="avg",logit="OR")
# summary(med.analysis,opt="avg",logit="logOR")
# summary(med.analysis,opt="avg",logit="all")
}
\keyword{datasets}
