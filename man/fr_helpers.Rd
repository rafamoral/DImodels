\name{fr_helpers} 
\alias{fr_helpers}
\alias{plot.fr_tost}
\alias{print.fr_tost}
\alias{summary.fr_tost}
\alias{delta_sd}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Helper functions for testing functional redundancy
}
\description{

These are helpers for the \code{\link{fr_test}} and \code{\link{fr_test_DI}} functions.

}

%}

\usage{
\method{print}{fr_tost}(x, ...)
\method{plot}{fr_tost}(x, xlab = "Difference in predictions", ...)
\method{summary}{fr_tost}(object, max.print = 20, ...)

delta_sd(k)
}


%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The data.frame returned by the \code{\link{fr_test}} or \code{\link{fr_test_DI}} functions.}
  \item{xlab}{The x-axis label for the forest plot display pairwise community contrast results.}
  \item{object}{The data.frame returned by the \code{\link{fr_test}} or \code{\link{fr_test_DI}} functions.}
  \item{max.print}{The maximum number of pairwise community contrasts to display in summary. Defaults to 20.}
  \item{k}{The multiplier to use for scaling residual model standard deviation in \code{delta_sd()}.}
  \item{...}{Additional arguments for base R generics. Ignored here.}
}

%%%%%%%%%%%%%%%%%%%%
%\details{
%details
%}

\value{

\code{print()} and \code{summary()} return the same data.frame passed in \code{object}.

\code{plot()} returns \code{NULL}.

\code{delta_sd()} return an object with class \code{delta_sd}.
}
\author{
Rishabh Vishwakarma, Matthias Suter, James O’Malley, John A. Finn, Rafael A. Moral, and Caroline Brophy
}

\seealso{
\code{\link{fr_test}}
\code{\link{fr_test_DI}}
\code{\link{DI}}
\code{\link{autoDI}}
}

\examples{

## Load the Switzerland data
  data(Switzerland)
## Summarise the Switzerland data
  summary(Switzerland)

## Fit a DI model
  m1 <- DI(y = "yield", prop = 4:7, treat = 'nitrogen', DImodel = 'AV',
  density = 'density', estimate_theta = TRUE, data = Switzerland)
  summary(m1)
  
## Test functional redundancy between species p1 and p2 
## use wrapper function for DI models
  t1 <- fr_test_DI(object = m1, redundant = c("p1", "p2"),
                   delta = 2)

## The overall test result along with results for failing criteria are printed
  print(t1)
## We fail to establish functional redundancy between p1 and p2 because their
## monoculture performances and p1:p2 interaction were not found to be equivalent

## The summary() function can be used to get p-values for each individual contrast
  summary(t1)

## A plot() function is also provided to visualise the test results
  plot(t1)

## Using delta_sd
  delta_sd(0.25)

}