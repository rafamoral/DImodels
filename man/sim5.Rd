\name{sim5}
\alias{sim5}
\encoding{UTF-8}
\docType{data}
\title{
The Simulated "sim5" Dataset
} 
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
The \code{sim5} dataset was simulated. There are nine species that vary in proportions (\code{p1 - p9}). It is assumed that species 1 to 5 come from functional group 1, species 6 and 7 from functional group 2 and species 8 and 9 from functional group 3. The response was simulated assuming that there were species identity effects and functional group specific interaction effects, with theta  = 0.7. 
}

\usage{data(sim5)} % all info below to be edited
\format{
  A data frame with 206 observations on the following 12 variables:
  \describe{
    \item{\code{community}}{A numeric vector identifying each unique community, i.e., two rows with the same community value also share the same set of p1 to p9 values.}
     \item{\code{richness}}{A numeric vector identifying the number of species in the initial composition.}
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
 
The dataset \code{sim5} was simulated with:  
  \itemize{
   \item identity effects for the nine species with values = 10, 9, 8, 7, 11, 6, 5, 8, 9
   \item functional group specific interaction effects; assume functional groups are labelled FG1, FG2 and FG3, then the interaction parameter values are: between FG1 and FG2 = 8, between FG1 and FG3 = 3, between FG2 and FG3 = 6, within FG1 = 6, within FG2 = 4 and within FG3 = 5 
   \item theta = 0.7 (where \eqn{\theta} is a non-linear parameter included as a power on each \eqn{pipj} product within interaction variables, see Connolly et al 2013 for details)
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
## Code to simulate the sim5 dataset
  
\donttest{

## Simulate dataset sim5 with 9 species and three functional groups.
## The species 1-5 are FG1, species 6-7 are FG2 and species 8-9 are FG3. 
## Assume ID effects and the FG interactions model, with theta = 0.7.

## Set up proportions
  data("design_a")
  sim5 <- design_a

## Create the functional group interaction variables, with theta = 0.7.
  FG_matrix <- DI_data(prop = 3:11, FG = c("FG1","FG1","FG1","FG1","FG1","FG2","FG2","FG3","FG3"),
                        data = sim5, theta = 0.7, what = "FG")
  sim5 <- data.frame(sim5, FG_matrix)
  names(sim5)[12:17] <- paste0(names(sim5)[12:17], "_theta")

## To simulate the response, first create a matrix of predictors that includes p1-p9, the 
##  treatment and the interaction variables.
  X <- model.matrix(~ p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 
                    + bfg_FG1_FG2_theta + bfg_FG1_FG3_theta + bfg_FG2_FG3_theta 
                    + wfg_FG1_theta + wfg_FG2_theta + wfg_FG3_theta -1, data = sim5)

## Create a vector of 'known' parameter values for simulating the response.
## The first nine are the p1-p9 parameters, and the second set of six are the interaction 
##  parameters. 
  sim5_coeff <- c(10,9,8,7,11, 6,5, 8,9,     8,3,6, 6,4,5)

##Create response and add normally distributed error 
  sim5$response <- as.numeric(X \%*\% sim5_coeff)
  set.seed(35748)
  r <- rnorm(n = 206, mean = 0, sd = 1.2)
  sim5$response <- round(sim5$response + r, digits = 3)
  sim5[,12:17] <- NULL
  
}



###########################
## Analyse the sim5 dataset

## Load the sim5 data
  data(sim5)
## View the first few entries
  head(sim5)
## Explore the variables in sim5
  str(sim5)

## Check characteristics of sim5
  hist(sim5$response)
  summary(sim5$response)
  plot(sim5$richness, sim5$response)
  plot(sim5$p1, sim5$response)
  plot(sim5$p2, sim5$response)
  plot(sim5$p3, sim5$response)
  plot(sim5$p4, sim5$response)
  plot(sim5$p5, sim5$response)
  plot(sim5$p6, sim5$response)
  plot(sim5$p7, sim5$response)
  plot(sim5$p8, sim5$response)
  plot(sim5$p9, sim5$response)

## What model fits best? Selection using F-test in autoDI 
  auto1 <- autoDI(y = "response", prop = 3:11, 
                  FG = c("FG1","FG1","FG1","FG1","FG1","FG2","FG2","FG3","FG3"), 
                  data = sim5, selection = "Ftest")
  summary(auto1)

## Fit the functional group model, with theta, using DI and the FG tag
  m1 <- DI(y = "response", prop = 3:11, 
           FG = c("FG1","FG1","FG1","FG1","FG1","FG2","FG2","FG3","FG3"), DImodel = "FG", 
           estimate_theta = TRUE, data = sim5)
  summary(m1)
  CI_95 <- theta_CI(m1, conf = .95)
  CI_95
  plot(m1)

\donttest{
## Check goodness-of-fit using a half-normal plot with a simulated envelope
  library(hnp)
  hnp(m1)
}

## Graph the profile likelihood
  library(ggplot2)
  ggplot(m1$profile_loglik, aes(x = grid, y = prof)) +
    theme_bw() +
    geom_line() +
    xlim(0,1.5) +
    xlab(expression(theta)) +
    ylab("Log-likelihood") + 
    geom_vline(xintercept = CI_95, lty = 3) + 
    labs(title = "   Log-likelihood versus theta", 
            caption = "dotted vertical lines are upper and lower bounds of 95\% CI for theta")
    
## Fit the functional group model, with theta set equal to the estimate from m1, and custom_formula.
## Note, it is not possible to estimate theta with custom_formula (only select a 'known' value). 
## First, create the functional group interactions (theta value as estimated from m1), 
##  store them in a new dataset and rename them with a theta indicator.
  FG_matrix <- DI_data(prop = 3:11, FG=c("FG1","FG1","FG1","FG1","FG1","FG2","FG2","FG3","FG3"), 
                      theta = 0.7296887, data = sim5, what = "FG")
  sim5new <- data.frame(sim5, FG_matrix)
  names(sim5new)[13:18] <- paste0(names(sim5new)[13:18], "_theta")
  m2 <- DI(custom_formula = response ~ 0 + p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 + 
         bfg_FG1_FG2_theta + bfg_FG1_FG3_theta + bfg_FG2_FG3_theta 
         + wfg_FG1_theta + wfg_FG2_theta + wfg_FG3_theta, data = sim5new)
## This will adjust the standard errors in m2 for the 'estimation' of theta
  m2$df.residual <- m2$df.residual - 1
## This will adjust the AIC in m2 for the 'estimation' of theta
  m2$aic <- m2$aic + 2
  summary(m2)

}

\keyword{datasets}
