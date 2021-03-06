% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blockmediate.R
\name{blockmediate}
\alias{blockmediate}
\title{blockmediate}
\usage{
blockmediate(
  lmodel.m,
  model.y,
  treat,
  treat.value = 1,
  control.value = 0,
  J = 1000,
  conf.level = 0.95,
  clust = NULL,
  seuil = 0
)
}
\arguments{
\item{lmodel.m}{list of fitted models object for mediators. Can be of class 'lm', 'polr','glm'.}

\item{model.y}{a fitted model object for the outcome. Can be of class 'lm', 'polr','glm'.}

\item{treat}{a character string indicating the name of the treatment variable used in the models. The treatment can be either binary (integer or a two-valued factor) or continuous (numeric).}

\item{treat.value}{value of the treatment variable used as the treatment condition. Default is 1.}

\item{control.value}{value of the treatment variable used as the control condition. Default is 0.}

\item{J}{number of Monte Carlo draws for quasi-Bayesian approximation.}

\item{conf.level}{level of the returned two-sided confidence intervals. Default is to return the 2.5 and 97.5 percentiles of the simulated quantities.}

\item{clust}{vector of length equal to the number of mediator. The order of the numbers corresponds to the order of the mediators in lmodel.m and the mediators belonging to the same cluster have the same number. Default is 'NULL', clusters are define using clustering spectral.}

\item{seuil}{threshold correlation value for spectral clustering. Default is 0.}
}
\value{
blockmediate returns an object of class "bm", a list that contains the components listed below.
The function summary (i.e., summary.bm) can be used to obtain a table of the results.
}
\description{
'blockmediate' is used to estimate various quantities for causal mediation analysis as "multimediate", not for several mediators but for for different block of mediators.
}
