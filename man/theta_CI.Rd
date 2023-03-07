\name{theta_CI}
\alias{theta_CI}
\encoding{UTF-8}
\title{
Compute Confidence Interval for Theta
}
\description{
This function allows the computation of a confidence interval for theta from a model object created from \code{DI} that includes the argument \code{estimate_theta = TRUE}, or from certain model objects created from \code{autoDI}.

A description of the non-linear parameter theta is available in Connolly et al 2013. 
}
\usage{
theta_CI(obj, conf = .95, n =100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{DI model object.}
  \item{conf}{Confidence level of the interval. The default is 0.95.}
  \item{n}{Number of subintervals in which the root is sought. The default is 100.}
}

\details{ 
The confidence interval calculated here is based on the values obtained when profiling the log-likelihood function for different values of theta. It is obtained in four steps:

1. define a grid of values for theta ranging from 0.01 to 2.5 of length 101 including the profile likelihood estimate of theta

2. fit the DI model setting theta equal to each value in the grid and obtain the log-likelihood value corresponding to each value of theta

3. obtain linear interpolations between the log-likelihood (\eqn{l}) values (here we use \code{\link{approxfun}})

4. calculate the lower and upper values of the CI by obtaining the values of theta corresponding to a log-likelihood value of \eqn{max_\theta(l) - 0.5*\chi^2_(1-\alpha;1)}, where \eqn{max_\theta(l)} is the maximum value of the profile log-likelihood obtained in the grid and \eqn{\chi^2_(conf;1)} is the \eqn{conf*100}\% percentile of the chi-squared distribution with 1 d.f.

When fitting any DI model setting \code{estimate_theta = TRUE}, steps 1 and 2 are automatically done within the \code{DI} function call. The \code{theta_CI} function performs steps 3 and 4 above to return the CI.

Note that when maximising the log-likelihood to find the estimate for theta, the parametric space is limited between 0.01 and 1.5. The larger grid (up to 2.5) is constructed to allow for obtaining the upper bound of the confidence interval in case the estimate of theta is close to 1.5.
}
\value{
The function returns a named numeric vector with two values: the \code{lower} and \code{upper} limits of the \code{conf}*100\% CI for theta.
}
\references{
Connolly J, T Bell, T Bolger, C Brophy, T Carnus, JA Finn, L Kirwan, F Isbell, J Levine, A \enc{Lüscher}{}, V Picasso, C Roscher, MT Sebastia, M Suter and A Weigelt (2013) An improved model to predict the effects of changing biodiversity levels on ecosystem function. Journal of Ecology, 101, 344-355.
}
\author{
Rafael A. Moral, John Connolly and Caroline Brophy
} 
%\note{
%further notes
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{DI}}
\code{\link{autoDI}}

Other examples using the \code{theta_CI} function:\cr
The \code{\link{Bell}} dataset examples.\cr
The \code{\link{sim2}} dataset examples.\cr
The \code{\link{sim5}} dataset examples.\cr
}
\examples{

## Load the sim5 data
  data(sim5)
## View the first five entries
  head(sim5)
## Explore the variables in sim5
  str(sim5)
  
## Fit the functional group model, with theta, using DI and the FG tag
  m1 <- DI(y = "response", prop = 3:11, 
           FG = c("FG1","FG1","FG1","FG1","FG1","FG2","FG2","FG3","FG3"), DImodel = "FG", 
           estimate_theta = TRUE, data = sim5)
  summary(m1)
  CI_95 <- theta_CI(m1, conf = .95)
  CI_95
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


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
