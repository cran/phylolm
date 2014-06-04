\name{phyloglm-methods}
\alias{print.phyloglm}
\alias{summary.phyloglm}
\alias{print.summary.phyloglm}
\alias{residuals.phyloglm}
\alias{vcov.phyloglm}
\alias{logLik.phyloglm}
\alias{print.logLik.phyloglm}
\alias{plot.phyloglm}
\title{Methods for class 'phyloglm'.}
\description{These are method functions for class 'phyloglm'.}
\usage{
\method{print}{phyloglm}(x, digits = max(3, getOption("digits") - 3), ...)
\method{summary}{phyloglm}(object, ...)
\method{residuals}{phyloglm}(object,type = c("response"), ...)
\method{vcov}{phyloglm}(object, ...)
\method{logLik}{phyloglm}(object, ...)
\method{plot}{phyloglm}(x, ...)
}
\arguments{
  \item{x}{an object of class \code{"phyloglm"}.}
  \item{object}{an object of class \code{"phyloglm"}.}
  \item{digits}{number of digits to show in summary method.}  
  \item{type}{Currently, only the "response" type is implemented. It returns the raw residuals, that is, the differences between the observed responses and the predicted values. They are phylogenetically correlated.}
  \item{\dots}{further arguments to methods.}
}

\author{Lam Si Tung Ho}
\seealso{
\code{\link{phyloglm}}
}
