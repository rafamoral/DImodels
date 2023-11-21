\name{predict.DI}
\alias{predict.DI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predict Method for Diversity-Interactions (DI) Models
}
\description{
Generates predictions for a fitted DI models object and, optionally, the associated standard errors for those predictions.
}
\usage{
\method{predict}{DI}(object, newdata, se.fit = FALSE,
interval = c("none", "confidence", "prediction"), 
level = 0.95, weights = 1,
type = c("link", "response", "terms"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a \code{\link{DI}} or \code{\link{autoDI}} model object.}
  \item{newdata}{optionally, a data frame in which to look for variables with which to predict. If omitted, predictions will be made for data used to fit the model.}
  \item{se.fit}{A logical switch indicating whether to calculate the associated standard errors.}
  \item{interval}{The type of interval to be calculated. Accepts one of "none", "confidence" or "prediction", with "none"" being the default.}
  \item{level}{The level (1 - \eqn{\alpha}) at which to calculate the interval. Defaults to 0.95, giving a 95\% interval.}
  \item{weights}{The variance weights for calculating the prediction intervals. By default all values get the same weight of 1. Can also specify a numeric vector of same length as number of rows in \code{newdata}.}
  \item{type}{the type of prediction required. One of "link", "response" or "terms". The default "link" is on the scale of the linear predictors while "response" is on the scale of the response variable. The "terms" option returns a matrix giving the fitted values of each term in the model formula on the linear predictor scale.}
  \item{...}{further arguments passed to or from other methods. For eg: \code{dispersion} or \code{na.action} arguments from \code{\link[stats]{predict.glm}} function.}
}
\details{
If \code{newdata} doesn't contain all predictors from the model, necessary numeric predictors will be added in with a value of 0, while categorical predictors will be added in with their baseline value. See examples for more information.
}
\value{
if \code{se.fit = FALSE}, a vector of predictions is returned.

if \code{se.fit = TRUE}, a list with the following components is returned.
  \item{fit}{Predictions, as for \code{se.fit = FALSE.}}
  \item{se.fit}{Estimated standard error for each prediction}
  \item{residual.scale}{A scalar giving the square root of the dispersion used in computing the standard errors.}
}
\references{
Connolly J, T Bell, T Bolger, C Brophy, T Carnus, JA Finn, L Kirwan, F Isbell, J Levine, A Luescher, V Picasso, C Roscher, MT Sebastia, M Suter and A Weigelt (2013) An improved model to predict the effects of changing biodiversity levels on ecosystem function. Journal of Ecology, 101, 344-355.

Kirwan L, J Connolly, JA Finn, C Brophy, A Luescher, D Nyfeler and MT Sebastia (2009) Diversity-interaction modelling - estimating contributions of species identities and interactions to ecosystem function. Ecology, 90, 2032-2038.
}
\author{
Rafael A. Moral, John Connolly, Rishabh Vishwakarma and Caroline Brophy
}


\seealso{
\code{\link{DI}}
\code{\link{autoDI}}
\code{\link[stats]{predict.glm}}
}

\examples{
## Load the Switzerland data
  data(Switzerland)
## Summarise the Switzerland data
  summary(Switzerland)

## Fit a DI model
  m1 <- DI(y = "yield", prop = 4:7, treat = 'nitrogen', DImodel = 'AV',
  density = 'density', estimate_theta = TRUE, data = Switzerland)
  summary(m1)
  
## Prediction without newdata, full dataset from model fit will be used
  predict(m1, se.fit = TRUE)

## Prediction with newdata
  newdata1 <- data.frame('p1' = c(1,0), 'p2' = c(0,1), 'p3' = c(0,0), 'p4' = c(0,0),
                         'nitrogen' = c(50, 150), 'density' = c('low','high'))
  predict(m1, newdata = newdata1, se.fit = TRUE)

## If any categorical variable is missing, the baseline level
## of the variable is used for prediction
  newdata2 <- newdata1[, -5]
  print(newdata2)
  predict(m1, newdata = newdata2)
  
## If a numerical variable is missing, the value 0 is used for it
  newdata3 <- newdata1[, -c(3,4)]
  print(newdata3)
  predict(m1, newdata = newdata3)
}
