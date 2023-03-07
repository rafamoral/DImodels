\name{contrasts_DI}
\alias{contrasts_DI}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Compute contrasts for Diversity Interactions (DI) models
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function calculates and tests contrasts using \code{\link[multcomp]{glht}} for a model object created by the \code{\link{DI}} or \code{\link{autoDI}} functions.
}
\usage{
contrasts_DI(object, contrast, contrast_vars, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A \code{\link{DI}} or \code{\link{autoDI}} model object.}
  \item{contrast}{The coefficients to generate the contrast. There are three options to specify these coefficients:
  
                  A list with each element being a vector of the same length as the number of linear coefficients in the model. The list can be named for tracking the contrasts.
                  
                  A matrix with the same number of columns as the number of linear coefficients in the model. The matrix can be given row names to track the different contrasts.
                  
                  A numeric vector with the length being a multiple of the number of coefficients in the matrix.
                  
                  Overrides \code{contrast_vars} if specified.}
  \item{contrast_vars}{A quicker and easier way to specify contrast coefficients for testing categorical variables. A nested named list with same names as the categorical variables in the model and same lengths as levels of those variables can be specified and the contrast will be calculated without the user having to specify the remaining variables. 
  
  Will be overridden if \code{contrast} is specified.
  
  Recommended for only testing cateogrical variables in the model.}
  
  \item{...}{Additional arguments passed to the \code{\link[multcomp]{glht}} function}
}
\details{
The contrasts are calculated and tested using the \code{\link[multcomp]{glht}} function in the \code{multcomp} package. 
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An object of class \code{\link[multcomp]{glht}} is returned
}
\references{
%% ~put references to the literature/web site here ~
Kirwan L, J Connolly, JA Finn, C Brophy, A \enc{Lüscher}{}, D Nyfeler and MT Sebastia (2009) Diversity-interaction modelling - estimating contributions of species identities and interactions to ecosystem function. Ecology, 90, 2032-2038. 

Torsten Hothorn, Frank Bretz and Peter Westfall (2008).
  Simultaneous Inference in General Parametric Models.
  Biometrical Journal 50(3), 346--363.
}
\author{
Rafael A. Moral, John Connolly, Rishabh Vishwakarma and Caroline Brophy
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{DI}}
\code{\link{autoDI}}
\code{\link[multcomp]{glht}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

## Load the Switzerland data
  data(Switzerland)
## Summarise the Switzerland data
  summary(Switzerland)

## Fit a DI model
  m1 <- DI(y = "yield", prop = 4:7, treat = 'nitrogen', DImodel = 'AV',
  density = 'density', estimate_theta = FALSE, data = Switzerland)
  summary(m1)
  
## Contrasts for difference between monocultures of p1 and p2, p3 and p4, p2 and p3
  con1 <- contrasts_DI(object = m1,
                       contrast = list('p1vp2 Mono' = c(1, -1, 0, 0, 0, 0, 0),
                                       'p3vp4 Mono' = c(0, 0, 1, -1, 0, 0, 0),
                                       'p2vp3 Mono' = c(0, -1, 1, 0, 0, 0, 0)))
  summary(con1)
  
## Contrasts for 50:50 mixture of p1 and p2 vs 50:50 mixture of p3 and p4  
  con2 <- contrasts_DI(object = m1,
                       contrast = list('p1p2 vs p3p4' = c(0.5, 0.5, -0.5, -0.5, 0, 0, 0)))
  summary(con2)

## Example using contrast_vars
  data(sim2)

### Fit model with block
  m2 <- DI(y = "response", prop = 3:6, DImodel = 'FULL', block = 'block',
           estimate_theta = FALSE, data = sim2)
  summary(m2)
  
### contrast for average of first two blocks vs third block
  con3 <- contrasts_DI(object = m2,
                       contrast_vars = list('block' = list('1_2vs3' = c(0.5, 0.5, -1, 0))))
  summary(con3)
}
