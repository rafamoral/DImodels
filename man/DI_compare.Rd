\name{DI_compare}
\alias{DI_compare}
\encoding{UTF-8}
\title{
Compare a Fitted Diversity-Interactions Model to the Reference Model
}
\description{
This function fits the reference model internally and compares a DI model fit to it using \code{anova}.
}
\usage{
DI_compare(model, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{A DI model object.}
  \item{...}{Further arguments passed to \code{anova}.}
}

\details{ 
This function takes a DI model object as input, fits the reference model internally and compares the two models using \code{anova}. The reference model includes an additive effect of community (each unique combination of species proportions) in the linear predictor. For more details on the reference model, see Connolly et al. (2013).
}
\value{
The function returns the reference model, a \code{glm} model object.
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
}
\examples{

## Load the sim1 data
  data(sim1)

## Fit the FULL model
  m1 <- DI(y = "response", block = "block", prop = 3:6, 
           DImodel = "FULL", data = sim1)

## Compare with the reference model
  DI_compare(m1, test = "F")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
