\name{autoDI}
\alias{autoDI}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Automated Diversity-Interactions Model Fitting
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function provides an automated way to fit a (limited) range of Diversity-Interactions (DI) models. Using one of several selection criteria, \code{autoDI} will identify the best DI model from the range fitted via a three-step selection process (see Details for more information). While \code{autoDI} can be a useful starting point for fitting DI models, its range of models is not exhaustive and additional model fitting or testing via \code{DI} is likely to be required. For instance, \code{autoDI} does not test for interactions of a treatment with other variables in the model.

}
\usage{
autoDI(y, prop, data, block, density, treat, ID, FG = NULL, 
       selection = c("Ftest","AIC","AICc","BIC","BICc"), 
       step0 = FALSE, step4 = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
%%     ~~Describe \code{y} here~~
The column name of the response vector, which must be in quotes, for example, \code{y = "yield"}. 
}
  \item{block}{
%%     ~~Describe \code{block} here~~
The name of the block variable (if present), which must be in quotes, for example, \code{block = "block"}. If no blocking variable, omit this argument.
}
  \item{density}{
%%     ~~Describe \code{density} here~~
The name of the density variable (if present), which must be in quotes, for example, \code{density = "density"}. If no density variable, omit this argument. 
}
  \item{prop}{
%%     ~~Describe \code{prop} here~~
A vector of \emph{s} column names identifying the species proportions in each row in the dataset. For example, if the species proportions columns are labelled p1 to p4, then \code{prop = c("p1","p2","p3","p4")}. Alternatively, the column numbers can be specified, for example, \code{prop = 4:7}, where the species proportions are in the 4th to 7th columns.
}
  \item{treat}{
%%     ~~Describe \code{prop} here~~
The name of a column in the dataset containing the value of a treatment factor or covariate. The treatment name must be included in quotes, for example, \code{treat = "nitrogen"}. (Only one treatment or covariate is permitted in \code{autoDI}, but see \code{DI} for options involving more than one treatment.) If the treatment is a factor, the variable must already be specified as a factor prior to using \code{autoDI}. 
}
  \item{ID}{
  %%     ~~Describe \code{ID} here~~
This argument takes a text list (of length \emph{s}) dsecirbing groupings for the identity effects of the species. For example, if there are four species and you wish to group the identity effects all four species into a single term: \code{ID} could be \code{ID = c("ID1","ID1","ID1","ID1")}, where "ID1" is the name of the ID group. Similarly if the we wish to have two identity effect groups where identity effect of species 1 and 3, and species 2 and 4 are grouped together: \code{ID} could be \code{ID = c("ID1","ID2","ID1","ID2")}, where "ID1" and "ID2" are the names of the ID groups. These ideas expand to any number of species and any number or combination of groups. Finally, the ID groups do not have to be named "ID1" and "ID2", the user can specify any name for the groups.

  \itemize{
   \item If the \code{ID} argument is not specified, each species will be assumed to have a separate identity effect. 
   
   \item Specify an grouping for the ID does not affect the interaction terms. The interactions are still calculated using the individual species proportions.

  \item The \code{ID} argument is defunct when using the \code{custom_formula} argument, since species identity effects must be included directly in the \code{custom_formula} argument. 
}
}
  \item{FG}{
%%     ~~Describe \code{FG} here~~
If species are classified by \emph{g} functional groups, this parameter gives a text list (of length \emph{s}) of the functional group to which each species belongs. For example, for four grassland species with two grasses and two legumes: FG could be \code{FG = c("G","G","L","L")}, where G stands for grass and L stands for legume.
}
  \item{data}{
%%     ~~Describe \code{data} here~~
Specify the dataset, for example, \code{data = Switzerland}. The dataset name should not appear in quotes.
}
\item{selection}{
%%     ~~Describe \code{selection} here~~
Selection method to be used in the automated model selection process. Options are \code{"Ftest"}, \code{"AIC"}, \code{"AICc"}, \code{"BIC"} and \code{"BICc"}. The default is \code{selection = "Ftest"}.
}
  \item{step0}{
%%     ~~Describe \code{step0} here~~
By default, \code{autoDI} outputs steps 1 - 4, however, an initial step 0 can be included by specifying \code{step0 = TRUE}. This will fit a model with an intercept only, and will sequentially add in and test the inclusion of \code{block}, \code{density} and \code{treat}, if they are specified in the \code{autoDI} arguments. 
}
  \item{step4}{
%%     ~~Describe \code{step4} here~~
Step 4 performs a lack of fit test for the final model selected by steps 1 - 3. By default, \code{step4 = TRUE}, but it can be omitted by specifying  \code{step4 = FALSE}.
}
}

%%%%%%%%%%%%%%%%%%%%
\details{
%%  ~~ If necessary, more details than the description above ~~

 
\strong{What are Diversity-Interactions models?} 
 

Diversity-Interactions (DI) models (Kirwan et al 2009) are a set of tools for analysing and interpreting data from experiments that explore the effects of species diversity on community-level responses. We recommend that users of the \code{DImodels} package read the short introduction to DI models (available at: \code{\link{DImodels}}). Further information on DI models is available in Kirwan et al 2009 and Connolly et al 2013.

 
\strong{Checks on data prior to using \code{autoDI}.} 
 

Before using \code{autoDI}, check that the species proportions for each row in your dataset sum to one. See the 'Examples' section for code to do this. An error message will be generated if the proportions don't sum to one.  

 
\strong{How does the \code{autoDI} function work?} 
 

The \code{autoDI} function identifies the 'best' Diversity-Interactions (DI) model from a specific range of proposed models using a three-step process of selection (Steps 1 to 3) and performs a lack of fit test on the selected model (Step 4). Only a limited subset of all possible models are examined under autoDI, for example, interactions involving experiment structural terms (\code{block}, \code{density}, \code{treat}) are not explored. The \code{autoDI} function can provide an excellent initial analysis, but often additional modelling and exploration will be required using the \code{DI} function. 

Steps 1-3 outlined below provide details on the automated model selection process followed by \code{autoDI}. Step 4 is a lack of fit test for the selected model. Step 0 may also be included (\code{step0 = TRUE}) as an initial step to sequentially test the inclusion of experiment structural variables (\code{block}, \code{density}, \code{treat}), prior to fitting the DI models in Step 1.

The default selection method used is \code{selection = "Ftest"}, which will return the appropriate F test statistic value(s) and p-value(s) in each step. When any of the information theoretic approaches (\code{"AIC"}, \code{"AICc"}, \code{"BIC"} or \code{"BICc"}) are specified, the model with the lowest value is selected in each step, even if it is only the lowest by a tiny margin; therefore, it is recommended to examine the information theoretic values across all models.

\emph{Step 1}

The AV model is fitted twice: considering the non-linear parameter \code{theta} equal to 1, and estimating it by maximising the profile likelihood. For example, for a design with 4 species, for the AV model, the parameter \code{theta} enters the model as:
 
\code{AV_theta = (p1*p2)^theta + (p1*p3)^theta + (p1*p4)^theta + (p2*p3)^theta + (p2*p4)^theta + (p3*p4)^theta}

In the profile likelihood estimation, \code{theta} is tested across the range 0.01 to 1.5. Then, the two models (with \code{theta} estimated and with \code{theta} set to 1) are compared using the method specified by the \code{selection} argument.

Details on DI models that include \code{theta} are described in Connolly et al 2013.


\emph{Step 2}

Five models are fitted, sequentially, each building on the previous model, and compared. If in Step 1 the conclusion is that \code{theta} is not different from 1, all models are fitted without estimating \code{theta}. If, however, in Step 1 the conclusion is that \code{theta} is different from 1, all models are fitted by fixing \code{theta} to its estimate obtained in Step 1.

If the experiment structural variables \code{treat}, \code{block} or \code{density} are specified, they will be included as an additive factor or covariate in \strong{each} of the five models, but interactions between them and the DI model terms will not be included or tested. 

Assume that \code{FG} (functional groups) have been specified. The five DI models are:

\itemize{
    \item \code{STR}: This model contains an intercept, and block, density and treatment if present; STR stand for 'structural' and represents experiment structural variables. 
    \item \code{ID}: The species proportions are added to the STR model. The proportions sum to 1, and are included in the model as \code{0 + p1 + ... + ps}, where \emph{s} is the number of species in the pool, as specified in the \code{prop} option.
    \item \code{AV}: The terms in the ID model, plus a single 'average' pairwise interaction term. For the \code{Switzerland} dataset, the single variable that is added to the model is computed as: \code{AV = p1*p2 + p1*p3 + p1*p4 + p2*p3 + p2*p4 + p3*p4}.
    \item \code{FG}: The terms in the ID model, plus interaction terms related to functional groups. These terms describe the average effects of interaction between pairs of species within each functional group, and between pairs of species from different functional groups. For example, in the \code{Switzerland} dataset there are  four species with \code{FG = c("Grass", "Grass", "Legume", "Legume")}, and there are six pairwise interactions, one between the two grasses, one between the two legumes, and four between grass and legume species. Grouping these interactions gives three terms, one for interactions between grasses, one for interactions between legumes and one for between a  grass and a legume species. The model assumes that any grass will interaction with any legume in the same way and the 'between functional group grass-legume' variable is computed as: \code{bfg_G_L = p1*p3 + p1*p4 + p2*p3 + p2*p4}. If there were more than two grasses in the dataset, this model would assume that any pair of grasses interact in the same way, similarly for legumes. If the \code{FG} argument is not specified, this model is omitted from Step 1.
    \item \code{FULL}: The terms in the ID model, plus an interaction term for each pair of species. When there are \emph{s} species, there are \emph{s*(s-1)/2} pairwise interaction terms, i.e., for the \code{Switzerland} dataset with four species, there are six interactions that are each added to the model (in R formula syntax: \code{p1:p2}, \code{p1:p3}, \code{p1:p4}, \code{p2:p3}, \code{p2:p4} and \code{p3:p4}).
}

If the \code{FG} argument is omitted, the \code{FG} model will be replaced by:
\itemize{
    \item \code{ADD}: The terms in the ID model, plus a species specific 'additive' interaction term for each species. These terms measure the interactive contribution of each species with any other species and are denoted \eqn{\lambda_i} for the \emph{i}th species. The interaction between any pair of species \emph{i} and \emph{j} is computed as \eqn{\lambda_i + \lambda_j}. 
}

If selection by information criteria is specified, both ADD and FG models will be fitted and compared using the selected information criterion.
 
\emph{Step 3}
 
If \code{treat} is specified, the selected model in Step 2 includes the \code{treat} variable (since all models in Step 2 include \code{treat} if present). Here, the selected model from Step 2 is re-fitted without \code{treat} and the models with and without the treatment are compared using the method specified by the \code{selection} argument.

If \code{treat} was not specified, this step is redundant.


 
 
\emph{Step 4} 

Step 4 provides a lack of fit test for the model selected by Steps 1 to 3. A factor called 'community' is created that has a level for each unique setting of the species proportions (as specified in the \code{prop} argument). The 'reference' model includes all terms in the model that was selected by Steps 1 to 3, plus the factor community. The reference model is compared to the model selected by Steps 1 to 3 via an F-test (regardless of the \code{selection} argument value), thus providing a lack of fit test. 

Note, that the reference model is never intended to be a candidate model, it is only fitted for the purpose of testing lack of fit.  If the test result is significant, it indicates that there is a lack of fit in the Diversity-Interactions model selected by \code{autoDI}.

Note also, that the test will not work if all combinations of the species proportions are unique. In this instance, the option \code{Step4 = FALSE} should be selected.   

}
\value{
The function returns the final selected model, an object of class \code{DI}.
}
\references{
%% ~put references to the literature/web site here ~
Connolly J, T Bell, T Bolger, C Brophy, T Carnus, JA Finn, L Kirwan, F Isbell, J Levine, A \enc{Lüscher}{}, V Picasso, C Roscher, MT Sebastia, M Suter and A Weigelt (2013) An improved model to predict the effects of changing biodiversity levels on ecosystem function. Journal of Ecology, 101, 344-355. 
  
Kirwan L, J Connolly, JA Finn, C Brophy, A \enc{Lüscher}{}, D Nyfeler and MT Sebastia (2009) Diversity-interaction modelling - estimating contributions of species identities and interactions to ecosystem function. Ecology, 90, 2032-2038. 
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

Other examples using \code{autoDI}: 
The \code{\link{sim1}} dataset examples. 
The \code{\link{sim2}} dataset examples. 
The \code{\link{sim3}} dataset examples. 
The \code{\link{sim4}} dataset examples. 
The \code{\link{sim5}} dataset examples. 
The \code{\link{Switzerland}} dataset examples. }
\examples{

## Load the Switzerland data
  data(Switzerland)
## Summarise the Switzerland data
  summary(Switzerland)
  
## Check that the proportions sum to 1 (required for DI models)
## p1 to p4 are in the 4th to 7th columns in Switzerland
  Switzerlandsums <- rowSums(Switzerland[4:7])
  summary(Switzerlandsums)
  
## Perform automated model fitting on the Switzerland dataset
  
## Model selection by F-test
  auto1 <- autoDI(y = "yield", density = "density", prop = c("p1","p2","p3","p4"), 
                  treat = "nitrogen", FG = c("G", "G", "L", "L"), data = Switzerland, 
                  selection = "Ftest")
  summary(auto1)
  
## Running autoDI after grouping identity effects using the "ID" argument
  auto2 <- autoDI(y = "yield", density = "density", 
                  prop = c("p1","p2","p3","p4"), 
                  treat = "nitrogen", ID = c("ID1", "ID1", "ID1", "ID1"),
                  FG = c("G", "G", "L", "L"), data = Switzerland, 
                  selection = "Ftest")
  summary(auto2)

  
\donttest{
## Using column numbers to indicate which columns contain the proportions and including Step 0
  auto2 <- autoDI(y = "yield", density = "density", prop = 4:7, treat = "nitrogen", 
                  FG = c("G", "G", "L", "L"), data = Switzerland, selection = "Ftest", step0 = TRUE)
  summary(auto2)
    
## Exclude the FG (functional group) argument to include the additive species "ADD" model in Step 1
  auto3 <- autoDI(y = "yield", density = "density", prop = 4:7, treat = "nitrogen", 
                  data = Switzerland, selection = "Ftest")
  summary(auto3)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
