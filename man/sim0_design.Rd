\name{sim0_design}
\alias{sim0_design}
\encoding{UTF-8}
\docType{data}
\title{
The "sim0_design" Dataset
} 
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
This dataset contains a set of proportions \code{p1} to \code{p3} where each row sums to 1. It is used as the design matrix for simulating other datasets in the \code{DImodels} package.  
}
\usage{data("sim0_design")}
\format{
  A data frame with 16 observations on the following 5 variables:
  \describe{
  \item{\code{community}}{A numeric vector identifying each unique community, i.e., two rows with the same community value also share the same set of p1 to p3 values.}
  \item{\code{richness}}{A numeric vector indicating the number of species in the initial composition, i.e., the number of proportions from p1 to p3 that are >0.}
    \item{\code{p1}}{A numeric vector indicating a proportion (of species 1).}
    \item{\code{p2}}{A numeric vector indicating a proportion (of species 2).}
    \item{\code{p3}}{A numeric vector indicating a proportion (of species 3).}
  }
}


\details{
%%  ~~ If necessary, more details than the __description__ above ~~
The columns p1 to p3 form a simplex space.
}


%\references{
%references
%}
\examples{

## Load the sim0_design data
  data(sim0_design)
## View the first five entries
  head(sim0_design)
## Explore the variables in design_a
  str(sim0_design)

## Check that the proportions sum to 1 (required for DI models)
## p1 to p3 are in the 3rd to 5th columns in sim0_design
  sim0_design_sums <- rowSums(sim0_design[3:5])
  summary(sim0_design)
  
}
\keyword{datasets}
