\name{DImodels-package}
\alias{DImodels-package}
\alias{DImodels}
\encoding{UTF-8}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\docType{package}
\title{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\packageTitle{DImodels}
}
\description{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\packageDescription{DImodels}
}
\details{

\strong{Introduction to Diversity-Interactions models}

Diversity-Interactions (DI) models (Kirwan et al 2009) are a set of tools for analysing and interpreting data from experiments that explore the effects of species diversity on community-level responses; for example, the effect of increasing community diversity on biomass production in a grassland ecosystem. Most analyses of diversity experiments quantify community diversity in terms of species richness, the number of species present. The DI method modifies this presence/absence approach in mixed communities by taking species relative abundance in the community into account. So, instead of ignoring differences in community responses across communities with the same species but with different species relative abundances, the DI approach aims to understand and explain these differences.

\strong{What variables will a suitable dataset contain?}

The DI approach assesses the effect of species diversity on a community response over a period of time. For data from a biodiversity study with \emph{S} species to be suited to the DI approach, the variables that are required on each experimental unit are:

\enumerate{
    \item A set of proportions \eqn{p1, p2, ... , pS} that characterise the proportions of each species at a defined starting point in time. These proportions (or relative abundances) of species in the community (\eqn{pi} for the ith species) range between 0 (absence) and 1 (monoculture - the only species present) and the sum of all the \eqn{pi} values for a community is always 1.   
    \item A community-level response variable, recorded a period of time after the initial species proportions were recorded.
  }
The dataset may also contain other variables such as a block, density or treatments.

\strong{What are Diversity-Interactions models?}

A DI model typically has three components and takes the form:

\deqn{y = Identities + Interactions + Structures + \epsilon}

where \emph{y} is a community-level response, the \emph{Identities} are the effects of species identities and enter the model as individual species proportions at the beginning of the time period, the \emph{Interactions} are the interactions among the species proportions, while \emph{Structures} include other experimental structures such as blocks, treatments or density. In a three species system, with an experimental blocking structure, a possible DI model is:

\deqn{y = \beta_1p1 + \beta_2p2 + \beta_3p3 + \delta_{12}p1p2 + \delta_{13}p1p3 + \delta_{23}p2p3 + \alpha_k + \epsilon}

Where \emph{p1}, \emph{p2} and \emph{p3} are the species proportions of species 1, 2 and 3 respectively, and \eqn{\alpha_k} is the effect of block \emph{k}. The error term \eqn{\epsilon} is usually assumed to be normally distributed with mean 0 and constant variance (and assumed independent and identically distributed). 

In a monoculture of the ith species, the expected performance of \emph{y} in block \emph{k} is \eqn{\beta_i}. In mixture, the \emph{Identities} component provides the expected performance as a weighted average of monoculture responses, and the \emph{Interactions} component is added to it to give the overall expected performance of the mixed community. The \emph{Interactions} component addresses the question: do mixed communities perform differently from what might be expected from the weighted averaging of monoculture performances? Note that in Kirwan et al 2009 the \emph{Interactions} component is referred to as the diversity effect, however, here we use the more general term "\emph{Interactions}", since the interpretation of how species diversity affects the response is a combination of both the \emph{Identities} and \emph{Interactions} components. The community with the best overall performance may depend on both the \emph{Identities} and \emph{Interactions} and will rarely be the community with the largest net interactions effect. The equation above provides an example of a DI model where there is a unique interaction term for each pair of species. It is possible to test various constraints among interactions, some of which may be motivated by the context of the data (Kirwan et al 2009). 

The \emph{Interactions} component may also include an non-linear exponent parameter \eqn{\theta} on each \eqn{pipj} term, where a value different to one allows the importance and impact of interaction terms to be altered (Connolly et al 2013). 

\strong{What can the \code{DImodels} package do?}

The DI approach is a full regression method where the response of the community is characterised by the effects of diversity variables such as  \emph{Identities} and \emph{Interactions} components, and by experimental structural variables such as blocks, density and treatments. All of these may be important determinants of response and so should be included in the analysis of community responses. The \code{DImodels} package aims to make it easier to analyse data using DI models. 

Currently, the \code{DImodels} package contains three main functions: \code{autoDI}, \code{DI} and \code{DI_data}. Here we give a brief overview of each, and link to the respective help files for further information. 

\code{auto_DI}: This function gives a simple overview of the successive contribution of the \emph{Structures}, \emph{Identities} and \emph{Interactions} components to the model via an automated model selection process. It will identify the 'best' model from a (limited) subset of all possible DI models. However, \code{autoDI} may need to be supplemented by more refined analysis, for example, \code{auto_DI} does not test for interactions between the terms in the \emph{Structures} and \emph{Identities} components. It can also only facilitate one block, one density and one treatment variable. However, it is a very useful starting point for DI model exploration. Further information at: \code{\link{autoDI}}.  

\code{DI}: This function can fit a wide range of DI models and includes the flexibility to test for multiple treatments or additional interactions, for example, between terms in the \emph{Identities} and \emph{Structures} components. The DI function fits one user-defined DI model at a time, and it allows for flexibility in the form of the model through a combination of in-built argument options and additional user-defined options. Further information at: \code{\link{DI}}. 

\code{DI_data}: This function can compute various forms of interactions among the \eqn{pi} variables. This function is not required when using \code{autoDI}, or for the in-built argument options in \code{DI}, but may be needed when specifying additional user-defined options in \code{DI}. Further information at: \code{\link{DI_data}} and for examples of when it is required see \code{\link{DI}}.  

There are three other functions:
\enumerate{
    \item \code{\link{theta_CI}} can fit a confidence interval to the parameter theta, when it has been estimated using either \code{autoDI} or \code{DI},
    \item \code{\link{DI_compare}} that can compare a fitted DI model to the 'reference' model (see \code{\link{autoDI}}) for details about the reference model), and
    \item \code{\link{update_DI}} that can update a fitted DI model passing different values to one or more arguments (see \code{\link{update_DI}}) for details).
}

\strong{Challenges with fitting and interpreting Diversity-Interactions models}

Analysing data using DI methods can be tricky for people familiar with ordinary regression models and the \code{DImodels} package aims to make analysis easier. The difficulties lie in:
\enumerate{
    \item The lack of familiarity in dealing with the \emph{Identities} component variables \eqn{pi}, which must sum to 1. This constraint can lead to interpretative issues, and estimation problems with some widely used R software.
    \item The novelty of specifying the \emph{Interactions} components in terms of many pairwise interaction terms whose numbers may greatly increase with \emph{S}. The number of pairwise terms can often be reduced by identifying biologically meaningful patterns among them, for example, through functional grouping of species. This may greatly reduce the number of coefficients to be estimated and interpreted in the model (Kirwan et al, 2009). 
    \item The introduction of a power coefficient theta (\eqn{\theta}) for all pairwise interaction terms (Connolly et al, 2013). This parameter can be very useful in describing the effect of community evenness on community response, i.e., whether response changes rapidly or slowly across communities where the relative abundance of species changes from being equal for all species to dominance by one or more species. 
}

\strong{Limitations of the \code{DImodels} package}

Currently, the \code{DImodels} package does not support multivariate responses or repeated measurements on the same experimental unit. It also does not currently support the random effects approach to modelling pairwise interactions that was developed by Brophy et al 2017.
}

\author{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\packageAuthor{DImodels}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Maintainer: \packageMaintainer{DImodels}
}
\references{
Brophy C, A Dooley, L Kirwan, JA Finn, J McDonnell, T Bell, MW Cadotte and J Connolly. (2017) Biodiversity and ecosystem function: Making sense of numerous species interactions in multi-species communities. Ecology 98, 1771-1778.

Connolly J, T Bell, T Bolger, C Brophy, T Carnus, JA Finn, L Kirwan, F Isbell, J Levine, A \enc{Lüscher}{}, V Picasso, C Roscher, MT Sebastia, M Suter and A Weigelt (2013) An improved model to predict the effects of changing biodiversity levels on ecosystem function. Journal of Ecology, 101, 344-355.

Kirwan L, J Connolly, JA Finn, C Brophy, A \enc{Lüscher}{}, D Nyfeler and MT Sebastia (2009) Diversity-interaction modelling - estimating contributions of species identities and interactions to ecosystem function. Ecology, 90, 2032-2038.
}
% Optionally other standard keywords, one per line, from file 
% KEYWORDS in the R documentation directory 
\keyword{ package }
\seealso{
\code{\link{autoDI}}
\code{\link{DI}}
\code{\link{DI_data}}
\code{\link{theta_CI}}

Example datasets:
The \code{\link{Bell}} dataset.
The \code{\link{sim1}} dataset.
The \code{\link{sim2}} dataset.
The \code{\link{sim3}} dataset.
The \code{\link{sim4}} dataset.
The \code{\link{sim5}} dataset.
The \code{\link{Switzerland}} dataset.
}
\examples{
## Load the Switzerland data
  data(Switzerland)
  
## Check that the proportions sum to 1 (required for DI models)
## p1 to p4 are in the 4th to 7th columns in Switzerland
  Switzerlandsums <- rowSums(Switzerland[4:7])
  summary(Switzerlandsums)
  
## Example of autoDI
  auto1 <- autoDI(y = "yield", density = "density", prop = c("p1","p2","p3","p4"), 
                  treat = "nitrogen", FG = c("G","G","L","L"), data = Switzerland, 
                  selection = "Ftest")
  summary(auto1)
  
## Example of DI
  m1 <- DI(y = "yield", density = "density", prop = 4:7, treat = "nitrogen", 
           FG = c("G","G","L","L"), DImodel = "FG", data = Switzerland)
  summary(m1)
  
## Example of DI_data.
## Create interaction variables and incorporate them into a new data frame Switzerland2.
## Switzerland2 will contain the new variables:  AV, E, p1_add, p2_add, p3_add, p4_add, 
##  bfg_G_L, wfg_G, wfg_L.
  FG_matrix <- DI_data(prop = c("p1","p2","p3","p4"), FG = c("G","G","L","L"),  
                       data = Switzerland, what = "FG")
  Switzerland2 <- data.frame(Switzerland, FG_matrix)
}
