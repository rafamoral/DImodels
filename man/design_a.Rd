\name{design_a}
\alias{design_a}
\encoding{UTF-8}
\docType{data}
\title{
The "design_a" Dataset
} 
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
This dataset contains a set of proportions \code{p1} to \code{p9} where each row sums to 1. It is used as the design matrix for simulating other datasets in the \code{DImodels} package.  
}
\usage{data("design_a")}
\format{
  A data frame with 206 observations on the following 11 variables:
  \describe{
  \item{\code{community}}{A numeric vector identifying each unique community, i.e., two rows with the same community value also share the same set of p1 to p9 values.}
  \item{\code{richness}}{A numeric vector indicating the number of species in the initial composition, i.e., the number of proportions from p1 to p9 that are >0.}
    \item{\code{p1}}{A numeric vector indicating a proportion (of species 1).}
    \item{\code{p2}}{A numeric vector indicating a proportion (of species 2).}
    \item{\code{p3}}{A numeric vector indicating a proportion (of species 3).}
    \item{\code{p4}}{A numeric vector indicating a proportion (of species 4).}
    \item{\code{p5}}{A numeric vector indicating a proportion (of species 5).}
    \item{\code{p6}}{A numeric vector indicating a proportion (of species 6).}
    \item{\code{p7}}{A numeric vector indicating a proportion (of species 7).}
    \item{\code{p8}}{A numeric vector indicating a proportion (of species 8).}
    \item{\code{p9}}{A numeric vector indicating a proportion (of species 9).}
  }
}


\details{
%%  ~~ If necessary, more details than the __description__ above ~~
The columns p1 to p9 form a simplex space.
}


%\references{
%references
%}
\examples{

## Load the design_a data
  data(design_a)
## View the first five entries
  head(design_a)
## Explore the variables in design_a
  str(design_a)

## Check that the proportions sum to 1 (required for DI models)
## p1 to p9 are in the 3rd to 11th columns in design_a
  design_a_sums <- rowSums(design_a[3:11])
  summary(design_a_sums)
  
}
\keyword{datasets}
