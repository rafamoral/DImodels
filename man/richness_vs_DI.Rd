\name{richness_vs_DI}
\alias{richness_vs_DI}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Comparing the Richness Model with DI Alternatives
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function provides an automated way to fit the richness model and a (limited) range of Diversity-Interactions (DI) models.

}
\usage{
richness_vs_DI(y, prop, data, extra_formula)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
%%     ~~Describe \code{y} here~~
The column name of the response vector, which must be in quotes, for example, \code{y = "yield"}. 
}
  \item{prop}{
%%     ~~Describe \code{prop} here~~
A vector of \emph{s} column names identifying the species proportions in each row in the dataset. For example, if the species proportions columns are labelled p1 to p4, then \code{prop = c("p1","p2","p3","p4")}. Alternatively, the column numbers can be specified, for example, \code{prop = 4:7}, where the species proportions are in the 4th to 7th columns.
}
  \item{data}{
%%     ~~Describe \code{data} here~~
Specify the dataset, for example, \code{data = Switzerland}. The dataset name should not appear in quotes.
}
\item{extra_formula}{
%%     ~~Describe \code{extra_formula} here~~
Additional terms can be added using \code{extra_formula}. A \code{~} must be included before specifying the terms. For example, \code{extra_formula = ~ treatment:AV} will add a treatment effect. Any variable included directly in \code{extra_formula} must already be contained in the dataset (interaction variables can be created using the function \code{\link{DI_data}}, if required).  
}
}

%%%%%%%%%%%%%%%%%%%%
\details{
%%  ~~ If necessary, more details than the description above ~~

Connolly et al. (2013; Appendix 1) shows that there is an equivalence between DI models and different types of richness models (linear and nonlinear predictors using richness in different scales).

This function fits four models and compares them using AIC. The four models are:

1. The richness model

\deqn{y = Intercept + Slope * Richness + \epsilon;}

2. The average pairwise interactions (AV) DI model with common identity effects and theta equal to 0.5 (which is equivalent to model 1 when communities are all equi-proportional);

3. The average pairwise interactions (AV) DI model allowing for unique identity effects, but maintaining theta equal to 0.5;

4. The average pairwise interactions (AV) DI model allowing for unique identity effects, and estimating theta.

The function prints a table with AIC and associated degrees of freedom for each of the four models above, and returns the model with the smallest AIC.

}
\value{
The function returns the final selected model, an object of class \code{DI} or \code{lm}.
}
\references{
%% ~put references to the literature/web site here ~
Connolly J, T Bell, T Bolger, C Brophy, T Carnus, JA Finn, L Kirwan, F Isbell, J Levine, A \enc{Lüscher}{}, V Picasso, C Roscher, MT Sebastia, M Suter and A Weigelt (2013) An improved model to predict the effects of changing biodiversity levels on ecosystem function. Journal of Ecology, 101, 344-355. 
}
\author{
Rafael A. Moral, John Connolly, Rishabh Vishwakarma and Caroline Brophy
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

## Load the Switzerland data
  data(Switzerland)


## compare the richness model with DI alternatives
  richness_vs_DI(y = "yield", prop = 4:7, data = Switzerland)

## include the density effects in the linear predictors of the four models
  richness_vs_DI(y = "yield", prop = 4:7, data = Switzerland, extra_formula = ~ density)
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
