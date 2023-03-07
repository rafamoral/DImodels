\name{design_b}
\alias{design_b}
\encoding{UTF-8}
\docType{data}
\title{
The "design_b" Dataset
} 
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
This dataset contains a set of proportions \code{p1} to \code{p6} where each row sums to 1. It is used as the design matrix for simulating other datasets in the \code{DImodels} package.  
}
\usage{data("design_b")}
\format{
  A data frame with 47 observations on the following seven variables:
  \describe{
    \item{\code{richness}}{A numeric vector indicating the number of species in the initial composition, i.e., the number of proportions from p1 to p6 that are >0.}
    \item{\code{p1}}{A numeric vector indicating a proportion (of species 1).}
    \item{\code{p2}}{A numeric vector indicating a proportion (of species 2).}
    \item{\code{p3}}{A numeric vector indicating a proportion (of species 3).}
    \item{\code{p4}}{A numeric vector indicating a proportion (of species 4).}
    \item{\code{p5}}{A numeric vector indicating a proportion (of species 5).}
    \item{\code{p6}}{A numeric vector indicating a proportion (of species 6).}
  }
}


\details{
%%  ~~ If necessary, more details than the __description__ above ~~
The columns p1 to p6 form a simplex space.
}


%\references{
%references
%}
\examples{

## Load the design_b data
  data(design_b)
## View the first five entries
  head(design_b)
## Explore the variables in design_b
  str(design_b)
  
## Check that the proportions sum to 1 (required for DI models)
## p1 to p6 are in the 2nd to 7th columns in design_b
  design_b_sums <- rowSums(design_b[2:7])
  summary(design_b_sums)
  
}
\keyword{datasets}
