\name{update_DI} 
\alias{update_DI}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Update DI function call
}
\description{

Update the call to the \code{\link{DI}} function. This function allows users to update a previous call to the \code{DI} function by only updating the relevant arguments from \code{DI}, without the need to write out the full \code{DI} code again.   

}

%}

\usage{
update_DI(object, ...)
}


%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a \code{DI} model object}
  \item{...}{other arguments passed to the \code{\link{DI}} function}
 
}

%%%%%%%%%%%%%%%%%%%%
%\details{
%details
%}

\value{
A model object of class \code{"glm"} including the components detailed in \code{\link{glm}}, plus the following:
\item{DIcall}{the call of the \code{DI} function}
\item{DIinternalcall}{an internal call made within the DI model fitting process}
}
\references{
%% ~put references to the literature/web site here ~
Connolly J, T Bell, T Bolger, C Brophy, T Carnus, JA Finn, L Kirwan, F Isbell, J Levine, A \enc{Lüscher}{}, V Picasso, C Roscher, MT Sebastia, M Suter and A Weigelt (2013) An improved model to predict the effects of changing biodiversity levels on ecosystem function. Journal of Ecology, 101, 344-355. 
  
Kirwan L, J Connolly, JA Finn, C Brophy, A \enc{Lüscher}{}, D Nyfeler and MT Sebastia (2009) Diversity-interaction modelling - estimating contributions of species identities and interactions to ecosystem function. Ecology, 90, 2032-2038. 
}
\author{
Rafael A. Moral, John Connolly and Caroline Brophy
}

\seealso{
\code{\link{DI}}
}

\examples{

## Load the Switzerland data
  data(Switzerland)

## Fit the FG DImodel, with factors density and treatment, and with theta = 1
  m1 <- DI(y = "yield", density = "density", prop = 4:7, treat = "nitrogen", 
         FG = c("G","G","L","L"), DImodel = "FG", data = Switzerland)
  summary(m1)

## Fit the FG DImodel, with factors density and treatment, and theta estimated
  m2 <- update_DI(m1, estimate_theta = TRUE)
  summary(m2)

## Fit the FULL DImodel, with factors density and treatment, and theta estimated
  m3 <- update_DI(m1, DImodel = "FULL", estimate_theta = TRUE)
  summary(m3)
}