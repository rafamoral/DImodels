\name{sim0}
\alias{sim0}
\encoding{UTF-8}
\docType{data}
\title{
The Simulated "sim0" Dataset
} 
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
The \code{sim0} dataset was simulated. There are four replicates and three species that vary in proportions (\code{p1 - p3}). There are 16 unique sets of proportions identified by the variable \code{community}. The response was simulated assuming that there were species identity effects and separate pairwise interactions effects.
}

\usage{data(sim0)}
\format{
  A data frame with 64 observations on the following six variables:
  \describe{
    \item{\code{community}}{A numeric vector identifying each unique community, i.e., two rows with the same community value also share the same set of p1 to p3 values.}
    \item{\code{richness}}{A numeric vector indicating the number of species in the initial composition, i.e., the number of proportions from p1 to p3 that are >0.}
    \item{\code{p1}}{A numeric vector indicating the initial proportion of species 1.}
    \item{\code{p2}}{A numeric vector indicating the initial proportion of species 2.}
    \item{\code{p3}}{A numeric vector indicating the initial proportion of species 3.}
    \item{\code{response}}{A numeric vector giving the simulated response variable.} 
  }
}


\details{
%%  ~~ If necessary, more details than the __description__ above ~~

\strong{What are Diversity-Interactions (DI) models?}

Diversity-Interactions (DI) models (Kirwan et al 2009) are a set of tools for analysing and interpreting data from experiments that explore the effects of species diversity on community-level responses. We strongly recommend that users read the short introduction to Diversity-Interactions models (available at: \code{\link{DImodels}}). Further information on Diversity-Interactions models is also available in Kirwan et al 2009 and Connolly et al 2013.

\strong{Parameter values for the simulation}

DI models take the general form of:

\deqn{y = Identities + Interactions + Structures + \epsilon}

where \emph{y} is a community-level response, the \emph{Identities} are the effects of species identities and enter the model as individual species proportions at the beginning of the time period, the \emph{Interactions} are the interactions among the species proportions, while \emph{Structures} include other experimental structures such as blocks, treatments or density.

The dataset \code{sim0} was simulated with:
  \itemize{
   \item identity effects for the four species with values = 25, 20, 15
   \item all 3 pairwise interaction effects with values: 30, 20, 40 (for pairs of species 1-2, 1-3, and 2-3, respectively). 
   \item \eqn{\epsilon} assumed normally distributed with mean 0 and standard deviation 2.
              }

}

%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~

%}
\references{
Connolly J, T Bell, T Bolger, C Brophy, T Carnus, JA Finn, L Kirwan, F Isbell, J Levine, A \enc{Lüscher}{}, V Picasso, C Roscher, MT Sebastia, M Suter and A Weigelt (2013) An improved model to predict the effects of changing biodiversity levels on ecosystem function. Journal of Ecology, 101, 344-355.

Kirwan L, J Connolly, JA Finn, C Brophy, A \enc{Lüscher}{}, D Nyfeler and MT Sebastia (2009) Diversity-interaction modelling - estimating contributions of species identities and interactions to ecosystem function. Ecology, 90, 2032-2038.
}
\examples{

## load the sim0 dataset
  data(sim0)

## Find the best DI model using autoDI and F-test selection
  auto1 <- autoDI(y = "response", prop = c("p1","p2","p3"), data = sim0, 
                  selection = "Ftest")
  summary(auto1)
  
## Fit the FULL model using DI and the ID tag
  m1 <- DI(y = "response", prop = c("p1","p2","p3"), DImodel = "FULL", 
           data = sim0)
  summary(m1)
  plot(m1)

\donttest{
## Check goodness-of-fit using a half-normal plot with a simulated envelope
  library(hnp)
  hnp(m1)
}

}
\keyword{datasets}
