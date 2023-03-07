\name{sim2}
\alias{sim2}
\encoding{UTF-8}
\docType{data}
\title{
The Simulated "sim2" Dataset
} 
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
The \code{sim2} dataset was simulated. There are four blocks and four species that vary in proportions (\code{p1 - p4}). There are 15 unique sets of proportions identified by the variable \code{community}. Each unique community appears once in each block. The response was simulated assuming that there were species identity effects, block effects, an average pairwise interaction effect and a theta value of 0.5. 
}

\usage{data(sim2)}
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
 
The dataset \code{sim2} was simulated with:  
  \itemize{
   \item identity effects for the four species with values = 10, 9, 8, 7
   \item block effects for the four blocks with values = 1, 1.5, 2, 0  
   \item an average pairwise interaction effect = 8
   \item theta = 0.5 (where \eqn{\theta} is a non-linear parameter included as a power on each \eqn{pipj} product within interaction variables, see Connolly et al 2013 for details)
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
## Code to simulate the sim2 dataset
  
\donttest{

## Simulate dataset sim2 with species identity effects, block effects and 
##  an average pairwise interaction effect with theta=0.5.

## Use the proportions from the first fifteen plots in Switzerland
  data(Switzerland)

## Repeat the 15 plots over four blocks. 
## Give each community type a unique (community) number.
  sim2 <- data.frame(community = rep(1:15, each = 4),
                   block = factor(rep(1:4, times = 15)),
                   p1 = rep(Switzerland$p1[1:15], each = 4),
                   p2 = rep(Switzerland$p2[1:15], each = 4),
                   p3 = rep(Switzerland$p3[1:15], each = 4),
                   p4 = rep(Switzerland$p4[1:15], each = 4))

## Create the average pairwise interaction variable, with theta = 0.5
  AV_variable <- DI_data(prop = c("p1","p2","p3","p4"), data = sim2, 
                         theta = 0.5, what = "AV")
  sim2 <- data.frame(sim2, "AV_theta" = AV_variable)

## To simulate the response, first create a matrix of predictors that includes p1-p4 and  
##  the four block variables and the average pairwise interaction variable with theta=0.5.
  X <- model.matrix(~ p1 + p2 + p3 + p4 + block + AV_theta -1, data = sim2)

## Create a vector of 'known' parameter values for simulating the response.
## The first four are the p1-p4 parameters, the second four are the block effects and 
##  the last one is the interaction parameter.
  sim2_coeff <- c(10,9,8,7,   1,1.5,2,0,    8)

## Create response and add normally distributed error 
  sim2$response <- as.numeric(X \%*\% sim2_coeff)
  set.seed(328781)
  r <- rnorm(n = 60, mean = 0, sd = 1.1)
  sim2$response <- round(sim2$response + r, digits = 3)
  sim2$AV_theta <- NULL
}



###################################################################################################
###################################################################################################
## sim2

###########################
## Analyse the sim2 dataset
  
## Load the sim2 data
  data(sim2)
## View the first few entries
  head(sim2)
## Explore the variables in sim2
  str(sim2)
  
## Check that the proportions sum to 1 (required for DI models)
## p1 to p4 are in the 3rd to 6th columns in sim2
  sim2sums <- rowSums(sim2[3:6])
  summary(sim2sums)
  
## Check characteristics of sim2
  hist(sim2$response)
  summary(sim2$response)
  plot(sim2$p1, sim2$response)
  plot(sim2$p2, sim2$response)
  plot(sim2$p3, sim2$response)
  plot(sim2$p4, sim2$response)

## Find the best DI model using autoDI and F-test selection
  auto1 <- autoDI(y = "response", prop = c("p1", "p2", "p3", "p4"), block = "block", data = sim2, 
                  selection = "Ftest")
  summary(auto1)
  
  ## Fit the average pairwise model, including theta, using DI and the AV tag
  m1 <- DI(y = "response", prop = c("p1","p2","p3","p4"), block = "block", DImodel = "AV", 
           estimate_theta = TRUE, data = sim2)
  summary(m1)
  CI_95 <- theta_CI(m1, conf = .95)
  CI_95
  plot(m1)
  library(hnp)

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
         
## Fit the average pairwise model, including theta, using DI and custom_formula
## A value of theta must be 'chosen'. Take: 0.4533437 from m1. The 'estimate_theta' option is not
##  available with custom_formula.
  AV_variable <- DI_data(prop = c(3:6), data = sim2, theta = 0.4533437, what = "AV")
  sim2a <- data.frame(sim2, "AV_theta" = AV_variable)
  m2 <- DI(y = "response", custom_formula = response ~ 0 + p1 + p2 + p3 + p4 + AV_theta + block, 
           data = sim2a)
## This will adjust the standard errors in m2 for the 'estimation' of theta
  m2$df.residual <- m2$df.residual - 1
## This will adjust the AIC in m2 for the 'estimation' of theta
  m2$aic <- m2$aic + 2
  summary(m2)
  
}
\keyword{datasets}
