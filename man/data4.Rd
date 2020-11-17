\name{data4}
\alias{data4}
\docType{data}
\title{
Toy dataset 4
}
\description{A toy dataset for causal mediation analysis with 3 uncausally related mediators.
}
\usage{data("data4")}
\format{
  A data frame with 1000 observations on the following 8 variables.
  \describe{
    \item{\code{Treatment}}{a numeric vector}
    \item{\code{C1}}{a factor with levels \code{0} \code{1}}
    \item{\code{C2}}{a factor with levels \code{0} \code{1}}
    \item{\code{C3}}{a factor with levels \code{0} \code{1}}
    \item{\code{M1}}{a numeric vector}
    \item{\code{M2}}{a factor with levels \code{0} \code{1}}
    \item{\code{M3}}{a factor with levels \code{0} \code{1} \code{2} \code{3}}
    \item{\code{Outcome}}{a numeric vector}
  }
}

\examples{
#data(data4)
#data4$Treatment=as.factor(data4$Treatment)
#data4$C1=as.factor(data4$C1)
#data4$C2=as.factor(data4$C2)
#data4$C3=as.factor(data4$C3)
#data4$M1=as.numeric(data4$M1)
#data4$M3=as.factor(data4$M3)
#data4$M2=as.factor(data4$M2)
#data4$Outcome=as.numeric(data4$Outcome)
#summary(data4)

#M1reg=lm(M1~  Treatment + C1 + C2 + C3, data = data4)
#M2reg=glm(M2~ Treatment + C1 + C3, data = data4, family = binomial("probit"))
#M3reg=polr(M3~Treatment + C2 + C3     , data = data4, method = "probit")

#Yreg=lm(Outcome~ Treatment + M1 + M2 + M3 + C1 + C2 + C3, data=data4)

#med.analysis=multimediate(lmodel.m=list(M1reg,M2reg,M3reg),
#correlated=TRUE,model.y=Yreg,treat="Treatment",
#treat.value=1,control.value=0,J=1000,conf.level=0.95)

#summary(med.analysis,opt="avg")
}
\keyword{datasets}
