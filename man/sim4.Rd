\name{sim4}
\alias{sim4}
\encoding{UTF-8}
\docType{data}
\title{
The Simulated "sim4" Dataset
} 
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
The \code{sim4} dataset was simulated. There is a covariate treatment and six species that vary in proportions (\code{p1 - p6}). It is assumed that species 1 and 2 come from functional group 1, species 3 and 4 from functional group 2 and species 5 and 6 from functional group 3. The response was simulated assuming that there were species identity effects, separate pairwise interaction effects and a covariate effect. 
}

\usage{data(sim4)} % all info below to be edited
\format{
  A data frame with 141 observations on the following nine variables:
  \describe{
    \item{\code{richness}}{A numeric vector identifying the number of species in the initial composition.}
    \item{\code{treatment}}{A covariate taking values 50, 150 or 250.}
    \item{\code{p1}}{A numeric vector indicating the initial proportion of species 1.}
    \item{\code{p2}}{A numeric vector indicating the initial proportion of species 2.}
    \item{\code{p3}}{A numeric vector indicating the initial proportion of species 3.}
    \item{\code{p4}}{A numeric vector indicating the initial proportion of species 4.}
    \item{\code{p5}}{A numeric vector indicating the initial proportion of species 5.}
    \item{\code{p6}}{A numeric vector indicating the initial proportion of species 6.}
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
 
The dataset \code{sim4} was simulated with:  
  \itemize{
   \item identity effects for the six species with values = 25, 16, 18, 20, 10, 12
   \item a covariate effect = 0.03
   \item all 15 pairwise interaction effects with values: 30, 27, 20, 15, 10, 9, 14, 18, 36, 17, 26, 32, 9, 21, 16 (for pairs of species 1-2, 1-3, 1-4, 1-5, 1-6, 2-3, 2-4, ... , 5-6 respectively). 
   \item theta = 1 (where \eqn{\theta} is a non-linear parameter included as a power on each \eqn{pipj} product within interaction variables, see Connolly et al 2013 for details)
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

####################################
## Code to simulate the sim4 dataset
  
\donttest{

## Simulate dataset sim4 with 6 species, three functional groups and three levels of a covariate
## The species 1-2 are FG1, species 3-4 are FG2 and species 5-6 are FG3. 
## Assume ID effects and the full pairwise interaction model, with a covariate.

## Set up proportions
  data("design_b")
  sim4a <- design_b

# Replicate the design for three values of a covariate
  sim4b <- sim4a[rep(seq_len(nrow(sim4a)), times = 3), ]
  sim4c <- data.frame(treatment = rep(c(50, 150, 250), each = 47))
  sim4 <- data.frame(richness = sim4b[,1], sim4c, sim4b[,2:7])
  row.names(sim4) <- NULL

## To simulate the response, first create a matrix of predictors that includes p1-p6, the treatment 
##  and all pairwise interaction variables 
  X <- model.matrix(~ p1 + p2 + p3 + p4 + p5 + p6 + treatment + (p1 + p2 + p3 + p4 + p5 + p6)^2 -1,
                    data = sim4)

## Create a vector of 'known' parameter values for simulating the response.
## The first six are the p1-p6 parameters, and the second set of one is the treatment parameter 
##  and the third set of 15 are the interaction parameters.
  sim4_coeff <- c(25,16,18,20,10,12,    0.03,      30,27,20,15,10,9,14,18,36,17,26,32,9,21,16)

## Create response and add normally distributed error 
  sim4$response <- as.numeric(X \%*\% sim4_coeff)
  set.seed(34261)
  r <- rnorm(n = 141, mean = 0, sd = 2)
  sim4$response <- round(sim4$response + r, digits = 3)

}


###########################
## Analyse the sim4 dataset

## Load the sim4 data
  data(sim4)
## View the first few entries
  head(sim4)
## Explore the variables in sim4
  str(sim4)

## Check characteristics of sim4
  hist(sim4$response)
  summary(sim4$response)
  plot(sim4$richness, sim4$response)
  plot(sim4$richness[sim4$treatment==50], sim4$response[sim4$treatment==50], ylim=c(0,40))
  plot(sim4$richness[sim4$treatment==150], sim4$response[sim4$treatment==150], ylim=c(0,40))
  plot(sim4$richness[sim4$treatment==250], sim4$response[sim4$treatment==250], ylim=c(0,40))
  plot(sim4$p1, sim4$response)
  plot(sim4$p2, sim4$response)
  plot(sim4$p3, sim4$response)
  plot(sim4$p4, sim4$response)
  plot(sim4$p5, sim4$response)
  plot(sim4$p6, sim4$response)

\donttest{
## What model fits best? Selection using F-test  
  auto1 <- autoDI(y = "response", prop = 3:8, treat = "treatment", 
                  FG = c("FG1","FG1","FG2","FG2","FG3","FG3"), data = sim4, selection = "Ftest")
  summary(auto1)

## Ignore functional groups (will replace FG model with ADD model in Step 1 selection)  
  auto2 <- autoDI(y = "response", prop = 3:8, treat = "treatment", data = sim4, selection = "Ftest")
  summary(auto2)
}

## Fit the functional group model using DI and the FG tag
  m1 <- DI(y = "response", prop = 3:8, treat = "treatment", 
           FG = c("FG1","FG1","FG2","FG2","FG3","FG3"), DImodel = FG, data = sim4)
  summary(m1)

## Fit the additive species model using DI and the ADD tag
  m2 <- DI(y = "response", prop = 3:8, treat = "treatment", DImodel = ADD, data = sim4)
  summary(m2)

## Fit the full pairwise model using DI and the FULL tag
  m3 <- DI(y = "response", prop = 3:8, treat = "treatment", DImodel = FULL, data = sim4)
  summary(m3)
  plot(m3)
  
\donttest{
## Check goodness-of-fit using a half-normal plot with a simulated envelope
  library(hnp)
  hnp(m3)
}

## Create the functional group and additive species interaction variables,
##  and store in a new data frame called sim4a
  newlist <- DI_data(prop = 3:8, FG = c("FG1","FG1","FG2","FG2","FG3","FG3"), 
                     data = sim4, what = c("FG", "ADD"))
  sim4a <- data.frame(sim4, newlist$FG, newlist$ADD)

## Fit the functional group model using DI and custom_formula (equivalent to m1)
  m4 <- DI(custom_formula = response ~ 0 + p1 + p2 + p3 + p4 + p5 + p6 + bfg_FG1_FG2 
           + bfg_FG1_FG3 + bfg_FG2_FG3 + wfg_FG1 + wfg_FG2 + wfg_FG3 + treatment, data = sim4a)
  summary(m4)

## Fit the additive species model using DI and custom_formula (equivalent to m2)
  m5 <- DI(custom_formula = response ~ 0 + p1 + p2 + p3 + p4 + p5 + p6 + p1_add 
           + p2_add + p3_add + p4_add + p5_add + p6_add + treatment, data = sim4a)
  summary(m5)

## Fit the full pairwise model using DI and custom_formula (equivalent to m3)
  m6 <- DI(custom_formula = response ~ 0 + p1 + p2 + p3 + p4 + p5 + p6 
           + (p1 + p2 + p3 + p4 + p5 + p6)^2 + treatment, data = sim4a)
  summary(m6)
  
}
\keyword{datasets}
