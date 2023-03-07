\name{sim3}
\alias{sim3}
\encoding{UTF-8}
\docType{data}
\title{
The Simulated "sim3" Dataset
} 
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
The \code{sim3} dataset was simulated. There are two treatments and nine species that vary in proportions (\code{p1 - p9}). It is assumed that species 1 to 5 come from functional group 1, species 6 and 7 from functional group 2 and species 8 and 9 from functional group 3. The response was simulated assuming that there were species identity effects and functional group specific interaction effects.
}

\usage{data(sim3)} % all info below to be edited
\format{
  A data frame with 412 observations on the following 13 variables:
  \describe{
    \item{\code{community}}{A numeric vector identifying each unique community, i.e., two rows with the same community value also share the same set of p1 to p9 values.}
     \item{\code{richness}}{A numeric vector identifying the number of species in the initial composition.}
    \item{\code{treatment}}{A factor with levels \code{A} or  \code{B}.}
    \item{\code{p1}}{A numeric vector indicating the initial proportion of species 1.}
    \item{\code{p2}}{A numeric vector indicating the initial proportion of species 2.}
    \item{\code{p3}}{A numeric vector indicating the initial proportion of species 3.}
    \item{\code{p4}}{A numeric vector indicating the initial proportion of species 4.}
    \item{\code{p5}}{A numeric vector indicating the initial proportion of species 5.}
    \item{\code{p6}}{A numeric vector indicating the initial proportion of species 6.}
    \item{\code{p7}}{A numeric vector indicating the initial proportion of species 7.}
    \item{\code{p8}}{A numeric vector indicating the initial proportion of species 8.}
    \item{\code{p9}}{A numeric vector indicating the initial proportion of species 9.}
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
 
The dataset \code{sim3} was simulated with:  
  \itemize{
   \item identity effects for the nine species with values = 10, 9, 8, 7, 11, 6, 5, 8, 9
   \item treatment effects = 3, 0
   \item functional group specific interact effects; assume functional groups are labelled FG1, FG2 and FG3, then the interaction parameter values are: between FG1 and FG2 = 4, between FG1 and FG3 = 9, between FG2 and FG3 = 3, within FG1 = 2, within FG2 = 3 and within FG3 = 1 
   \item theta = 1 (where \eqn{\theta} is a non-linear parameter included as a power on each \eqn{pipj} product within interaction variables, see Connolly et al 2013 for details)
   \item \eqn{\epsilon} assumed normally distributed with mean 0 and standard deviation 1.2.
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
## Code to simulate the sim3 dataset
  
\donttest{

## Simulate dataset sim3 with 9 species, three functional groups and two levels of a treatment.
## The species 1-5 are FG1, species 6-7 are FG2 and species 8-9 are FG3. 
## Assume ID effects and the FG interaction model, with a treatment (factor with two levels).

## Set up proportions
  data("design_a")
  sim3a <- design_a

# Replicate the design over two treatments
  sim3b <- sim3a[rep(seq_len(nrow(sim3a)), each = 2), ]
  sim3c <- data.frame(treatment = factor(rep(c("A","B"), times = 206))) 
  sim3 <- data.frame(sim3b[,1:2], sim3c, sim3b[,3:11])
  row.names(sim3) <- NULL

## Create the functional group interaction variables
  FG_matrix <- DI_data(prop = 4:12, FG = c("FG1","FG1","FG1","FG1","FG1","FG2","FG2","FG3","FG3"),
                        data = sim3, what = "FG")
  sim3 <- data.frame(sim3, FG_matrix)

## To simulate the response, first create a matrix of predictors that includes p1-p9, the treatment 
##  and the interaction variables. 
  X <- model.matrix(~ p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 + treatment 
                    + bfg_FG1_FG2 + bfg_FG1_FG3 + bfg_FG2_FG3 
                    + wfg_FG1 + wfg_FG2 + wfg_FG3 -1, data=sim3)

## Create a vector of 'known' parameter values for simulating the response.
## The first nine are the p1-p9 parameters, and the second set of two are the treatment effects 
##   and the third set of six are the interaction parameters.
  sim3_coeff <- c(10,9,8,7,11, 6,5, 8,9,   3,0,      4,9,3, 2,3,1)

## Create response and add normally distributed error 
  sim3$response <- as.numeric(X \%*\% sim3_coeff)
  set.seed(1657914)
  r <- rnorm(n = 412, mean = 0, sd = 1.2)
  sim3$response <- round(sim3$response + r, digits = 3)

  sim3[,13:18] <- NULL
}



###########################
## Analyse the sim3 dataset

## Load the sim3 data
  data(sim3)
## View the first few entries
  head(sim3)
## Explore the variables in sim3
  str(sim3)

## Check characteristics of sim3
  hist(sim3$response)
  summary(sim3$response)
  plot(sim3$richness, sim3$response)
  plot(sim3$p1, sim3$response)
  plot(sim3$p2, sim3$response)
  plot(sim3$p3, sim3$response)
  plot(sim3$p4, sim3$response)
  plot(sim3$p5, sim3$response)
  plot(sim3$p6, sim3$response)
  plot(sim3$p7, sim3$response)
  plot(sim3$p8, sim3$response)
  plot(sim3$p9, sim3$response)

\donttest{
## What model fits best? Selection using F-test in autoDI  
  auto1 <- autoDI(y = "response", prop = 4:12, treat = "treatment", 
                  FG = c("FG1","FG1","FG1","FG1","FG1","FG2","FG2","FG3","FG3"), data = sim3, 
                  selection = "Ftest")
  summary(auto1)
}

## Fit the functional group model, with treatment, using DI and the FG tag
  m1 <- DI(y = "response", prop = 4:12, 
           FG = c("FG1","FG1","FG1","FG1","FG1","FG2","FG2","FG3","FG3"), treat = "treatment", 
           DImodel = "FG", data = sim3)
  summary(m1)
  plot(m1)

\donttest{
## Check goodness-of-fit using a half-normal plot with a simulated envelope
  library(hnp)
  hnp(m1)
}

## Create the functional group interactions and store them in a new dataset
  FG_matrix <- DI_data(prop = 4:12, FG = c("FG1","FG1","FG1","FG1","FG1","FG2","FG2","FG3","FG3"), 
                        data = sim3, what = "FG")
  sim3a <- data.frame(sim3, FG_matrix)

## Test if the FG interaction variables interact with treatment using 'extra_formula'
  m2 <- DI(y = "response", prop = 4:12, 
           FG = c("FG1","FG1","FG1","FG1","FG1","FG2","FG2","FG3","FG3"),
           treat = "treatment", DImodel = "FG", extra_formula = ~ bfg_FG1_FG2:treatment 
           + bfg_FG1_FG3:treatment + bfg_FG2_FG3:treatment + wfg_FG1:treatment + wfg_FG2:treatment 
           + wfg_FG3:treatment, data = sim3a)
  summary(m2)

## Fit the functional group model using DI and custom_formula
## Set up a dummy variable for treatment first (required).
  sim3a$treatmentA <- as.numeric(sim3a$treatment=="A")
  m3 <- DI(y = "response", custom_formula = response ~ 0 + p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 
           + p9 + treatmentA + bfg_FG1_FG2 + bfg_FG1_FG3 + bfg_FG2_FG3 + wfg_FG1 + wfg_FG2 
           + wfg_FG3, data = sim3a)
  summary(m3)
}
\keyword{datasets}
