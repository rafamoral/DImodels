\name{extract} 
\alias{extract}
\alias{extract.DI}
\alias{extract.autoDI}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract terms from a DI or autoDI model object
}
\description{

Extracts terms from a \code{\link{DI}} or \code{\link{autoDI}} model object that were calculated using the data preparation functions prior to the model fitting procedure.

}

%}

\usage{
extract(obj, what = c("FULL","ADD","FG","AV","E"))
\method{extract}{DI}(obj, what = c("FULL","ADD","FG","AV","E"))
\method{extract}{autoDI}(obj, what = c("FULL","ADD","FG","AV","E"))

}


%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{a \code{DI} or \code{autoDI} model object}
  \item{what}{the variable(s) to be extracted from the object.}
 
}

%%%%%%%%%%%%%%%%%%%%
%\details{
%details
%}

\value{
A list containing one or more of the following elements:
\item{E}{the evenness interaction variable}
\item{AV}{the average pairwise interaction variable}
\item{FG}{the functional group-related interaction variables}
\item{ADD}{the additive species contributions interaction variables}
\item{FULL}{all pairwise interactions}
}
\references{
%% ~put references to the literature/web site here ~
Connolly J, T Bell, T Bolger, C Brophy, T Carnus, JA Finn, L Kirwan, F Isbell, J Levine, A \enc{Lüscher}{}, V Picasso, C Roscher, MT Sebastia, M Suter and A Weigelt (2013) An improved model to predict the effects of changing biodiversity levels on ecosystem function. Journal of Ecology, 101, 344-355. 
  
Kirwan L, J Connolly, JA Finn, C Brophy, A \enc{Lüscher}{}, D Nyfeler and MT Sebastia (2009) Diversity-interaction modelling - estimating contributions of species identities and interactions to ecosystem function. Ecology, 90, 2032-2038. 
}
\author{
Rafael A. Moral, John Connolly and Caroline Brophy
}

\seealso{
\code{\link{DI}}
}

\examples{

## Load the Switzerland data
  data(Switzerland)

## Fit the FG DImodel, with factors density and treatment, and with theta = 1
  m1 <- DI(y = "yield", density = "density", prop = 4:7, treat = "nitrogen", 
           FG = c("G","G","L","L"), DImodel = "FG", data = Switzerland)

## Extract the FG terms
  extract(m1, what = "FG")

}