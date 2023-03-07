\name{sim1}
\alias{sim1}
\encoding{UTF-8}
\docType{data}
\title{
The Simulated "sim1" Dataset
} 
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
The \code{sim1} dataset was simulated. There are four blocks and four species that vary in proportions (\code{p1 - p4}). There are 15 unique sets of proportions identified by the variable \code{community}. Each unique community appears once in each block. The response was simulated assuming that there were species identity effects and block effects, but no diversity effects.
}

\usage{data(sim1)}
\format{
  A data frame with 60 observations on the following seven variables:
  \describe{
    \item{\code{community}}{A numeric vector identifying each unique community, i.e., two rows with the same community value also share the same set of p1 to p4 values.}
    \item{\code{block}}{A factor taking values 1 to 4 indicating block membership.}
    \item{\code{p1}}{A numeric vector indicating the initial proportion of species 1.}
    \item{\code{p2}}{A numeric vector indicating the initial proportion of species 2.}
    \item{\code{p3}}{A numeric vector indicating the initial proportion of species 3.}
    \item{\code{p4}}{A numeric vector indicating the initial proportion of species 4.}
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

The dataset \code{sim1} was simulated with:
  \itemize{
   \item identity effects for the four species with values = 10, 9, 8, 7
   \item block effects for the four blocks with values = 1, 1.5, 2, 0  
   \item no interaction effects
   \item \eqn{\epsilon} assumed normally distributed with mean 0 and standard deviation 1.1.
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

####################################
## Code to simulate the sim1 dataset

\donttest{

## Simulate dataset sim1 with species identity effects and block effects, but no interaction effect

## Use the proportions from the first fifteen plots in Switzerland
  data(Switzerland)

## Repeat the 15 plots over four blocks.
## Give each community type a unique (community) number.
  sim1 <- data.frame(community = rep(1:15, each = 4),
                   block = factor(rep(1:4, times = 15)),
                   p1 = rep(Switzerland$p1[1:15], each = 4),
                   p2 = rep(Switzerland$p2[1:15], each = 4),
                   p3 = rep(Switzerland$p3[1:15], each = 4),
                   p4 = rep(Switzerland$p4[1:15], each = 4))

## To simulate the response, first create a matrix of predictors that includes 
##  p1-p4 and the four block dummy variables.
  X <- model.matrix(~ p1 + p2 + p3 + p4 + block -1, data = sim1)

## Create a vector of 'known' parameter values for simulating the response.
## The first four are the p1-p4 parameters, the second four are the block effects.
  sim1_coeff <- c(10,9,8,7,   1,1.5,2,0)

## Create response and add normally distributed error 
  sim1$response <- as.numeric(X \%*\% sim1_coeff)
  set.seed(2020)
  r <- rnorm(n = 60, mean = 0, sd = 1.1)
  sim1$response <- round(sim1$response + r, digits = 3)
}



###########################
## Analyse the sim1 dataset

## Load the sim1 data
  data(sim1)
## View the first few entries
  head(sim1)
## Explore the variables in sim1
  str(sim1)

## Check that the proportions sum to 1 (required for DI models)
## p1 to p4 are in the 3rd to 6th columns in sim1
  sim1sums <- rowSums(sim1[3:6])
  summary(sim1sums)
  
## Check characteristics of sim1
  hist(sim1$response)
  summary(sim1$response)
  plot(sim1$p1, sim1$response)
  plot(sim1$p2, sim1$response)
  plot(sim1$p3, sim1$response)
  plot(sim1$p4, sim1$response)
  
## Find the best DI model using autoDI and F-test selection
  auto1 <- autoDI(y = "response", prop = c("p1","p2","p3","p4"), block = "block", data = sim1, 
                  selection = "Ftest")
  summary(auto1)
  
## Fit the identity model using DI and the ID tag
  m1 <- DI(y = "response", prop = c("p1","p2","p3","p4"), block = "block", DImodel = "ID", 
           data = sim1)
  summary(m1)
  plot(m1)

\donttest{
## Check goodness-of-fit using a half-normal plot with a simulated envelope
  library(hnp)
  hnp(m1)
}
  
## Fit the identity model using DI and custom_formula
  m2 <- DI(y = "response", custom_formula = response ~ 0 + p1 + p2 + p3 + p4 + block, data = sim1)
  summary(m2)

}
\keyword{datasets}
