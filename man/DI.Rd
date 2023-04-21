\name{DI} 
\alias{DI}
\alias{DI_ADD_treat}
\alias{DI_ADD}
\alias{DI_FULL_treat}
\alias{DI_FULL}
\alias{DI_AV_treat}
\alias{DI_AV}
\alias{DI_CUSTOM}
\alias{DI_E_treat}
\alias{DI_E}
\alias{DI_FG_treat}
\alias{DI_FG}
\alias{DI_ID_treat}
\alias{DI_ID}
\alias{DI_STR_treat}
\alias{DI_STR}
\alias{DI_theta}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Diversity-Interactions Model Fitting
}
\description{

This function will fit a wide range of Diversity-Interactions (DI) models, one at a time. It provides some assisted automated ways to fit DI models, and includes the flexibility to extend DI models in several directions.  

}

%}

\usage{
DI(y, prop, DImodel, custom_formula, data,
   block, density, treat, FG, extra_formula,
   estimate_theta = FALSE, theta = 1)
}


%- maybe also 'usage' for other objects documented here.
\arguments{
The minimum required arguments to use DI are either: 
\itemize{
    \item Argument \code{DImodel} with \code{data}, \code{y} and \code{prop}, or 
    \item Argument \code{custom_formula} with \code{data}.
  }
The \code{DImodel} argument allows fitting of DI models via a range of 'tag' options that determine the form of the species interactions terms (the tags, described below, are \code{STR}, \code{ID}, \code{AV}, \code{FG}, \code{ADD} and \code{FULL}) and extra terms can be added to the model using the \code{extra_formula} argument. Using the argument \code{custom_formula} requires full specification of the model to be fitted using standard \code{lm} or \code{glm} syntax. 
 

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
A vector of \emph{s} column names identifying the species proportions in each community in the dataset. For example, if the species proportions columns are labelled p1 to p4, then \code{prop = c("p1","p2","p3","p4")}. Alternatively, the column numbers can be specified, for example, \code{prop = 4:7}, where species proportions are in the 4th to 7th columns.
}
  \item{treat}{
%%     ~~Describe \code{treat} here~~
The name of a column in the dataset containing the value of a treatment factor or covariate. The treatment name must be included in quotes, for example, \code{treat = "nitrogen"}. If the treatment is a factor, the variable must already be specified as a factor prior to using \code{DI}. 

  \itemize{
   \item 
When used in conjunction with \code{DImodel}, the treatment will be included in the model as an additive factor or covariate, for example, specifying \code{treat = nitrogen, DImodel = ID} will fit the model \code{p1 + p2 + ... + ps + nitrogen}. Additional treatments, or interactions between the treatment and other model terms can be included via the \code{extra_formula} argument.

\item The \code{treat} argument is defunct when using the \code{custom_formula} argument, and any treatment must be included directly in the \code{custom_formula} argument.
}
}
  \item{FG}{
  %%     ~~Describe \code{FG} here~~
If species are classified by \emph{g} functional groups, this argument takes a text list (of length \emph{s}) of the functional group to which each species belongs. For example, for four grassland species with two grasses and two legumes: FG could be \code{FG = c("G","G","L","L")}, where G stands for grass and L stands for legume.

  \itemize{
   \item The \code{FG} argument is required if \code{DImodel = "FG"} is specified. 

  \item The \code{FG} argument is defunct when using the \code{custom_formula} argument, since species interactions must be included directly in the \code{custom_formula} argument. 
}
}
  \item{DImodel}{
%%     ~~Describe \code{DImodel} here~~
This argument is chosen (over \code{custom_formula}) to fit an automated version of a DI model. The chosen tag determines the form of the species interactions to be included in the model. The tags (or options) are: 
  \itemize{
   \item \code{STR} (no identity or interaction effects, only an intercept is fitted, plus the experiment structural variables \code{block}, \code{density} and \code{treat}, if specified).  
   }
   Each of the following includes the species proportions as specified in \code{prop}, the interaction variables according to the tag, plus \code{block}, \code{density} and \code{treat} if specified. 
   \itemize{
   \item \code{ID} (no interaction terms),  
   \item \code{AV} (a single average species pairwise interaction variable),  
   \item \code{FG} (functional group interaction variables, the \code{FG} argument must be specified to use this option),  
   \item \code{ADD} (the additive species interaction variables),  
   \item \code{FULL} (all pairwise interaction variables). 
   }
   The DImodel tag should appear in quotes, for example, \code{DImodel = "STR"}.
}
\item{extra_formula}{
%%     ~~Describe \code{extra_formula} here~~
In conjunction with \code{DImodel}, additional terms can be added using \code{extra_formula}. A \code{~} must be included before specifying the terms. For example, if \code{DImodel = "AV"} has been specified, adding \code{extra_formula = ~ I(AV**2)} will add a quadratic average pairwise interaction variable or \code{extra_formula = ~ treatment:AV} will add an interaction between the average pairwise species interaction variable and the treatment. Any variable included directly in \code{extra_formula} must already be contained in the dataset (interaction variables can be created using the function \code{\link{DI_data}}, if required).  
}
  \item{custom_formula}{
%%     ~~Describe \code{custom_formula} here~~
To specify your own DI model, write your own model formula using the \code{custom_formula} argument. The standard notation from \code{lm} and \code{glm} is used here, for example, \code{custom_formula = yield ~ 0 + p1:treatment + p2:treatment + p3:treatment + p4:treatment} will fit the DI model with the identity effects each crossed with treatment, where treatment is a variable in the dataset. The \code{custom_formula} argument is recommended when the \code{DImodel} and \code{extra_formula} arguments are not sufficient. Any variable included directly in \code{custom_formula} must be already contained in the dataset (interaction variables can be created using the function \code{\link{DI_data}}, if required).  
}
  \item{data}{
%%     ~~Describe \code{data} here~~
Specify the dataset, for example, \code{data = Switzerland}. The dataset name should not appear in quotes.
}
  \item{estimate_theta}{
%%     ~~Describe \code{estimate_theta} here~~
By default, theta (the power parameter on all \eqn{pi*pj} components of each interaction variable in the model) is set equal to one. Specify \code{estimate_theta = TRUE} to include the estimation of \eqn{\theta} in the specified model. The \code{estimate_theta} argument can only be used in conjunction with the \code{DImodel} argument; if the \code{custom_formula} is used, then theta estimation is not available and the default \code{estimate_theta = FALSE} must be used. 
}
  \item{theta}{
  Users may specify a value of theta different than 1 to fit the DI model. Note that if \code{estimate_theta = TRUE}, then theta will be estimated via maximum profile log-likelihood and the value specified for \code{theta} will be overridden.
  }
 
}

%%%%%%%%%%%%%%%%%%%%
\details{

 
\strong{What are Diversity-Interactions models?} 
 

Diversity-Interactions (DI) models (Kirwan et al 2009) are a set of tools for analysing and interpreting data from experiments that explore the effects of species diversity on community-level responses. We recommend that users of the \code{DImodels} package read the short introduction to DI models (available at: \code{\link{DImodels}}). Further information on DI models is available in Kirwan et al 2009 and Connolly et al 2013.

 
\strong{Checks on data prior to using \code{DI}.} 
 

Before using \code{DI}, check that the species proportions for each row in your dataset sum to one. See the 'Examples' section for code to do this. An error message will be generated if the proportions don't sum to one.  


 
\strong{How does the \code{DI} function work?} 
 
The \code{DI} function provides wide flexibility in the types of Diversity-Interactions (DI) models that can be fitted. There are two ways to fit models in \code{DI}: 1) using \code{DImodel}, possibly augmented by \code{extra_formula}, or 2) using \code{custom_formula}. Models are estimated using iteratively reweighted least squares, via the \code{\link{glm}} package, when the option \code{estimate_theta = FALSE}. 

Consider the following DI model, for example (in R formula syntax):  
\code{y ~ p1 + p2 + p3 + treatment + p1:p2 + p1:p3 + p2:p3 + p1:p2:treatment + p1:p3:treatment + p2:p3:treatment} 
 
This model can be fitted using \code{DImodel} and \code{extra_formula}:  
\code{DI(y = "y", prop = c("p1", "p2", "p3"), treat = "nitrogen", DImodel = "FULL", extra_formula = ~ p1:p2:treatment + p1:p3:treatment + p2:p3:treatment, data = datasetname)} 
 
or, by specifying all of the terms in the model using \code{custom_formula}: 
\code{DI(custom_formula = y ~ p1 + p2 + p3 + treatment + p1:p2 + p1:p3 + p2:p3 + p1:p2:treatment + p1:p3:treatment + p2:p3:treatment, data = datasetname)} 
 
We recommend to use \code{DImodel} where possible, to augment with \code{extra_formula} where required, and to only use \code{custom_formula} when \code{DImodel} plus \code{extra_formula} is insufficient.  


 
\strong{Including theta in DI models} 
 
A non-linear parameter \eqn{\theta} can be included in DI models as a power on all \eqn{pi*pj} components of each pairwise interaction variable. For example (in R formula syntax): 
\code{y ~ p1 + p2 + p3 + (p1:p2)^theta + (p1:p3)^theta + (p2:p3)^theta}  
for the full pairwise interaction model. Including \eqn{\theta} alters the contribution of the interaction term to the response (Connolly et al 2013).  

By default, the value of \eqn{\theta} is 1. By specifying \code{estimate_theta = TRUE} within \code{DI}, a value of \eqn{\theta} will be estimated using profile likelihood over the space \eqn{\theta} = 0.01 to 1.5. The option \code{estimate_theta = TRUE} can only be used with \code{DImodel}, it is not available when using \code{custom_formula}.

As a general guideline to testing if \eqn{\theta} is required, we recommend:

1) finding the best form of the species interaction terms assuming theta = 1, and then,

2) testing if theta differs from 1.

If no species interaction terms are needed, then there is no need to do any testing for theta.
}

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
\code{\link{autoDI}} 
\code{\link{theta_CI}}

Other examples using \code{DI}: 
The \code{\link{Bell}} dataset examples. 
The \code{\link{sim1}} dataset examples. 
The \code{\link{sim2}} dataset examples. 
The \code{\link{sim3}} dataset examples. 
The \code{\link{sim4}} dataset examples. 
The \code{\link{sim5}} dataset examples. 
The \code{\link{Switzerland}} dataset examples. 
}

\examples{ 

## Load the Switzerland data
data(Switzerland)

## Check that the proportions sum to 1 (required for DI models)
## p1 to p4 are in the 4th to 7th columns in Switzerland
Switzerlandsums <- rowSums(Switzerland[4:7])
summary(Switzerlandsums)

## Fit the FG DImodel, with factors density and treatment, and with theta = 1
m1 <- DI(y = "yield", density = "density", prop = 4:7, treat = "nitrogen", 
         FG = c("G","G","L","L"), DImodel = "FG", data = Switzerland)
summary(m1)

## Fit the FG DImodel, with factors density and treatment, and theta estimated
m2 <- DI(y = "yield", density = "density", prop = 4:7, treat = "nitrogen", 
         FG = c("G","G","L","L"), DImodel = "FG", data = Switzerland, estimate_theta = TRUE)
summary(m2)

## Test if the identity effects interact with nitrogen (and main effect of nitrogen excluded)
m3 <- DI(y = "yield", density = "density", prop = 4:7, FG = c("G", "G", "L", "L"), DImodel = "FG", 
         extra_formula = ~ (p1 + p2 + p3 + p4):nitrogen, data = Switzerland)
summary(m3)

## Fit the average pairwise model and check for a quadratic term using extra_formula.
## Need to create AV variable to be included in extra_formula and put in new dataset Switzerland2.
  AV_variable <- DI_data(prop = c("p1","p2","p3","p4"), data = Switzerland, what = "AV")
  Switzerland2 <- data.frame(Switzerland, "AV" = AV_variable)
  m4 <- DI(y = "yield", density = "density", prop = 4:7, DImodel = "AV", 
           extra_formula = ~ I(AV**2), data = Switzerland2)
  summary(m4)

## Using the custom_formula option to fit some, but not all, of the FG interactions.
## Fit the FG DImodel using custom_formula, with factors density and treatment, and theta = 1.
## Need to create functional group interaction variables for inclusion in custom_formulaand put
##  in new dataset Switzerland3.
  FG_matrix <- DI_data(prop = 4:7, FG = c("G","G","L","L"), data = Switzerland, what = "FG")
  Switzerland3 <- data.frame(Switzerland, FG_matrix)
  m5 <- DI(y = "yield", prop = c("p1","p2","p3","p4"), 
           custom_formula = yield ~ 0 + p1 + p2 + p3 + p4 + bfg_G_L + wfg_G 
           + density + nitrogen,
           data = Switzerland3)
  summary(m5)

}