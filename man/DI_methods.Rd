\name{DI_methods}
\alias{anova.DI}
\alias{AIC.DI}
\alias{BIC.DI}
\alias{AICc.DI}
\alias{BICc.DI}
\alias{logLik.DI}
\alias{fortify.DI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Methods for \code{DI} Objects
}
\description{
Different methods that can be used with objects of class \code{DI}.
}
\usage{
\method{anova}{DI}(object, ...)
\method{AIC}{DI}(object, ...)
\method{BIC}{DI}(object, ...)
\method{AICc}{DI}(obj)
\method{BICc}{DI}(obj)
\method{logLik}{DI}(object, ...)
\method{fortify}{DI}(model, data, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a \code{DI} model object}
  \item{obj}{a \code{DI} model object}
  \item{model}{a \code{DI} model object}
  \item{data}{data to which to add model fit statistics. Defaults to the model frame.}
\item{...}{further arguments passed to \code{anova}, \code{AIC}, \code{BIC}, \code{loglik}, or \code{fortify}}
}
\author{
Rafael A. Moral, John Connolly, Rishabh Vishwakarma and Caroline Brophy
}

\seealso{
\code{\link{autoDI}}
}
%\examples{
%## examples
%}