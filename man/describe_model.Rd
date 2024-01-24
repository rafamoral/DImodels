\name{describe_model}
\alias{describe_model}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Describe DI models
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function accepts a DImodel (i.e., regression models fit using the \code{\link{DI}} or \code{\link{autoDI}} functions) object and returns a short text summary describing the model. 
}
\usage{
describe_model(model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{A \code{\link{DI}} or \code{\link{autoDI}} regression model object.}
}

\value{
A short text describing the supplied DImodel object.
}
\references{
%% ~put references to the literature/web site here ~
Kirwan L, J Connolly, JA Finn, C Brophy, A \enc{Lüscher}{}, D Nyfeler and MT Sebastia (2009) Diversity-interaction modelling - estimating contributions of species identities and interactions to ecosystem function. Ecology, 90, 2032-2038. 
}
\author{
Rafael A. Moral, John Connolly, Rishabh Vishwakarma and Caroline Brophy
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{DI}}
\code{\link{autoDI}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

## Load the  data
   data(sim2)
## Fit model
   mod_FG <- DI(y = "response", FG = c("G", "G", "L", "L"), 
                prop = 3:6, data = sim2, DImodel = "FG")
## Describe model
   describe_model(mod_FG)


   mod_FULL <- DI(y = "response", estimate_theta = TRUE, 
                  prop = 3:6, data = sim2, DImodel = "FULL")
   describe_model(mod_FULL)

   mod_AV <- DI(y = "response", ID = c("ID1", "ID1", "ID2", "ID2"),
                estimate_theta = TRUE, 
                prop = 3:6, data = sim2, DImodel = "AV")
   describe_model(mod_AV)
}
