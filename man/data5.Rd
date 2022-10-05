\name{data5}
\alias{data5}
\docType{data}
\title{
Toy dataset 5
}
\description{A toy dataset for causal mediation analysis with 3 uncausally related mediators
and survival outcome.
}
\usage{data("data4")}
\format{
  A data frame with 3000 observations on the following 8 variables.
  \describe{
    \item{\code{Exposure}}{a numeric vector}
    \item{\code{M1}}{a numeric vector}
    \item{\code{M2}}{a numeric vector}
    \item{\code{M3}}{a numeric vector}
    \item{\code{event}}{event}
    \item{\code{surv_time}}{time before event}
  }
}

\examples{
#data(data5)
#modM1=lm(M1 ~ Exposure, data = data5)
#modM2=lm(M2 ~ Exposure, data = data5)
#modM3=lm(M3 ~ Exposure, data = data5)
#lmodel.m=list(modM1,modM2, modM3)
#model.y=aalen(Surv(surv_time, event) ~ const(Exposure) + const(M1) + const(M2)
#+ const(M3), data = data5, robust=T)
#multi.media=multimediate(lmodel.m,correlated=TRUE,model.y,treat='Exposure',
#treat.value=1,control.value=0,J=1000,conf.level=0.95,data=data5)
#summary(multi.media, opt='avg')
}
\keyword{datasets}
