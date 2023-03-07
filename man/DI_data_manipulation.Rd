\name{DI_data_manipulation}
\alias{DI_data_manipulation}
\alias{DI_data}
\alias{DI_data_ADD}
\alias{DI_data_E_AV}
\alias{DI_data_FG}
\alias{DI_data_prepare}
\alias{DI_data_fullpairwise}
\encoding{UTF-8}
\title{
Data manipulation function
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function can be used to compute additional variables for the various types of interactions among pairs of species proportions. These variables are defined following Kirwan et al 2007 and 2009, and Connolly et al 2013. 

The variables are denoted:  
 
\code{AV} and \code{E}: the average pairwise interaction variable (AV) and a scaled version of it (E), each a single variable. The variable AV is the sum of products of the  proportions of each pair of species in the mixture. The variable E is a scaled version of AV that ranges between 0 (for a monoculture community) to 1 for the equi-proportional mixture of all species in the pool.  

\code{FG}: the functional group interaction variables. There is a variable for (within) each functional group and one for (between) each pair of functional groups, i.e., if there are two functional groups, there will be three functional group interaction variables, while if there are three functional groups, there will be six functional group interaction variables.  
 
\code{ADD}: the additive species interaction variables, one for each species. 

\code{FULL}: all individual pairwise interactions. There will be \eqn{s*(s-1)/2} new variables created, where \emph{s} is the number of species in the pool.


By default, the interaction variables described above are created with \code{theta = 1}, but a different value of theta can also be specified (Connolly et al 2013). 
}
\usage{
DI_data(prop, FG, data, theta = 1, what = c("E","AV","FG","ADD","FULL"))
}


%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prop}{ 
%%     ~~Describe \code{P} here~~
A vector of column names identifying the species proportions in the dataset. For example, if the species proportions columns are labelled p1 to p4, then \code{prop = c("p1","p2","p3","p4")}. The column numbers in which the proportions are stored can also be referred to, for example,  \code{prop = 4:7} for the \code{Switzerland} data.
}
  \item{FG}{
%%     ~~Describe \code{FG} here~~
If species are classified by \emph{g} functional groups, this argument takes a text list (of length \emph{s}) of the functional group to which each species belongs. For example, for four grassland species with two grasses and two legumes, it could be \code{FG = c("G","G","L","L")}, where G stands for grass and L stands for legume. This argument is required when \code{"FG"} is included in the \code{what} argument. 
}
  \item{data}{
%%     ~~Describe \code{data} here~~
Specify the dataset, for example, \code{data = Switzerland}. The dataset name should not appear in quotes. 
}
\item{theta}{
%%     ~~Describe \code{data} here~~
Interaction variables will be computed with the theta power, equal to the value specified, on all \eqn{pi*pj} components of each interaction variable, with default value one. For example, with three species \code{AV = p1*p2 + p1*p3 + p2*p3} and if computed with \code{theta} = 0.5, this becomes \code{(p1*p2)^0.5 + (p1*p3)^0.5 + (p2*p3)^0.5}.
}
\item{what}{
%%     ~~Describe \code{y} here~~
The interactions to be computed. By default each set of variables will be computed: AV, E, FG (if \code{FG} argument specified), ADD and FULL. Individual sets can be selected. For example, \code{what = c("FULL")} or \code{what = c("AV", "FULL")}.  
}


}
\details{ 
%%  ~~ If necessary, more details than the description above ~~

 
\strong{What are Diversity-Interactions models?} 
 

Diversity-Interactions (DI) models (Kirwan et al 2009) are a set of tools for analysing and interpreting data from experiments that explore the effects of species diversity on community-level responses. We recommend that users of the \code{DImodels} package read the short introduction to DI models (available at: \code{\link{DImodels}}). Further information on DI models is available in Kirwan et al 2009 and Connolly et al 2013.

 
\strong{Checks on data prior to using the \code{DI_data} function.} 
 

Before using the \code{DI_data} function, check that the species proportions in each row of your dataset sum to one. See the 'Examples' section for code to do this. An error message will be generated if the proportions don't sum to one. 

 
\strong{When is the \code{DI_data} function needed?} 
 

It is not required to use the \code{DI_data} function if using the \code{autoDI} function, or the \code{DImodel} option in the \code{DI} function, as they will automatically create the species interaction variables needed. If using species interaction variables in the \code{extra_formula} or \code{custom_formula} options in \code{DI}, then it is required to have the variables already in the dataset and \code{DI_data} can do that.  

 
\strong{Short worked example to illustrate how the \code{DI_data} function works} 
 

The code to implement this example is provided in the 'Examples' section. 

Assume four species with initial proportions in two communities: (0.1, 0.2, 0.3, 0.4) and (0.25, 0.25, 0.25, 0.25), with \code{FG = c("G","G","L","L")}. 
 
For community 1: (0.1,0.2,0.3,0.4), assuming theta = 1, the data preparation functions will compute the following additional variables (details in Kirwan et al 2007 and 2009): 

AV = 0.1*0.2 + 0.1*0.3 + 0.1*0.4 + 0.2*0.3 + 0.2*0.4 + 0.3*0.4 = 0.35

E = \emph{(2s/(s-1))}*AV = 0.9333

p1_add = 0.1 * (1 - 0.1) = 0.09

p2_add = 0.2 * (1 - 0.2) = 0.16

p3_add = 0.3 * (1 - 0.3) = 0.21

p4_add = 0.4 * (1 - 0.4) = 0.24

bfg_G_L = 0.1*0.3 + 0.1*0.4 + 0.2*0.3 + 0.2*0.4 = 0.21

wfg_G = 0.1*0.2 = 0.02

wfg_L = 0.3*0.4 = 0.12

For community 1: (0.1,0.2,0.3,0.4), assuming theta = 0.5, the data preparation functions will compute the follow additional variables (details in Connolly et al 2013): 

AV = (0.1*0.2)^0.5 + (0.1*0.3)^0.5 + (0.1*0.4)^0.5 + (0.2*0.3)^0.5 + (0.2*0.4)^0.5 + (0.3*0.4)^0.5 = 1.3888
 
E =\emph{(2s/(s-1))}*AV = 3.7035

p1_add = 0.1^0.5 * (0.2^0.5 + 0.3^0.5 + 0.4^0.5) = 0.5146

p2_add = 0.2^0.5 * (0.1^0.5 + 0.3^0.5 + 0.4^0.5) = 0.6692

p3_add = 0.3^0.5 * (0.1^0.5 + 0.2^0.5 + 0.4^0.5) = 0.7646

p4_add = 0.4^0.5 * (0.1^0.5 + 0.2^0.5 + 0.3^0.5) = 0.8293

bfg_G_L = (0.1*0.3)^0.5 + (0.1*0.4)^0.5 + (0.2*0.3)^0.5 + (0.2*0.4)^0.5 = 0.9010

wfg_G = (0.1*0.2)^0.5 = 0.1414

wfg_L = (0.3*0.4)^0.5 = 0.3464  

When using the \code{DI_data} function to create interactions for theta values for a value different from 1, it is recommended to rename the new interaction variables to include \code{_theta}. 

The data manipulation values for community 2 (0.25, 0.25, 0.25,0.25) can be seen when the 'Examples' section code is run.
}
\value{
\code{DI_data} returns a named list with one or more the following components (depending on the specification of the \code{what} argument:
\item{E}{a numeric vector containing the evenness interaction variable}
\item{AV}{a numeric vector containing the average pairwise interaction variable}
\item{FG}{a matrix containing the functional group-related interaction variables}
\item{ADD}{a matrix containing the additive species contributions interaction variables}
\item{FULL}{a matrix containing all pairwise interactions}
}
\references{
Connolly J, T Bell, T Bolger, C Brophy, T Carnus, JA Finn, L Kirwan, F Isbell, J Levine, A \enc{Lüscher}{}, V Picasso, C Roscher, MT Sebastia, M Suter and A Weigelt (2013) An improved model to predict the effects of changing biodiversity levels on ecosystem function. Journal of Ecology, 101, 344-355. 
 
Kirwan L, A \enc{Lüscher}{}, MT Sebastia, JA Finn, RP Collins, C Porqueddu, A Helgadottir, OH Baadshaug, C Brophy, C Coran, S Dalmannsdottir, I Delgado, A Elgersma, M Fothergill, BE Frankow-Lindberg, P Golinski, P Grieu, AM Gustavsson, M \enc{Höglind}{}, O Huguenin-Elie, C Iliadis, M \enc{Jørgensen}{}, Z Kadziuliene, T Karyotis, T Lunnan, M Malengier, S Maltoni, V Meyer, D Nyfeler, P Nykanen-Kurki, J Parente, HJ Smit, U Thumm, & J Connolly (2007) Evenness drives consistent diversity effects in intensive grassland systems across 28 European sites. Journal of Ecology, 95, 530-539.  
 
Kirwan L, J Connolly, JA Finn, C Brophy, A \enc{Lüscher}{}, D Nyfeler and MT Sebastia (2009) Diversity-interaction modelling - estimating contributions of species identities and interactions to ecosystem function. Ecology, 90, 2032-2038. 
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

Other examples using the \code{DI_data} function: 
The \code{\link{sim2}} dataset examples. 
The \code{\link{sim3}} dataset examples. 
The \code{\link{sim4}} dataset examples. 
The \code{\link{sim5}} dataset examples. 
The \code{\link{Switzerland}} dataset examples. 
}
\examples{
################################
  
#### Data manipulation for the Switzerland dataset
  
## Load the Switzerland data
  data(Switzerland)
  
## Check that the proportions sum to 1 (required for DI models)
## p1 to p4 are in the 4th to 7th columns in Switzerland
  Switzerlandsums <- rowSums(Switzerland[4:7])
  summary(Switzerlandsums)
  
  
## Create FG interaction variables and incorporate them into a new data frame Switzerland2.
## Switzerland2 will contain the new variables:  bfg_G_L, wfg_G and wfg_L.
  FG_matrix <- DI_data(prop = c("p1","p2","p3","p4"), FG = c("G","G","L","L"), 
                       data = Switzerland, what = "FG")
  Switzerland2 <- data.frame(Switzerland, FG_matrix)
  
  
## Create FG interaction variables and incorporate them into a new data frame Switzerland3.
## Use theta = 0.5.
  FG_matrix <- DI_data(prop = c("p1","p2","p3","p4"), FG = c("G","G","L","L"), 
                       data = Switzerland, what = "FG", theta = 0.5)
  Switzerland3 <- data.frame(Switzerland, FG_matrix)
## Add "_theta" to the new interaction variables to differentiate from when theta = 1
  names(Switzerland3)[9:11] <- paste0(names(Switzerland3)[9:11], "_theta") 


#### All interactions can be added to a new dataset all together:

## Create all interactions and add them to a new data frame called Switzerland4
  newlist <- DI_data(prop = c("p1","p2","p3","p4"), FG = c("G","G","L","L"), data = Switzerland, 
                     what = c("E", "AV", "FG", "ADD", "FULL"))
  Switzerland4 <- data.frame(Switzerland, "E" = newlist$E, "AV" = newlist$AV, newlist$FG, 
                             newlist$ADD, newlist$FULL)

  
#### Or the various interactions can also be added to a new dataset individually:
  
## Create the average pairwise interaction and evenness variables
##  and store them in a new data frame called Switzerland5.
## Switzerland5 will contain the new variables: AV, E
  E_variable <- DI_data(prop = c("p1","p2","p3","p4"), data = Switzerland, what = "E")
  AV_variable <- DI_data(prop = c("p1","p2","p3","p4"), data = Switzerland, what = "AV")
  Switzerland5 <- data.frame(Switzerland, "AV" = AV_variable, "E" = E_variable)
  
## Create the functional group variables and add them to Switzerland5.
## In the FG names vector: G stands for grass, L stands for legume.
## Switzerland5 will contain: bfg_G_L, wfg_G and wfg_L
  FG_matrix <- DI_data(prop = 4:7, FG = c("G","G","L","L"), data = Switzerland, what = "FG")
  Switzerland5 <- data.frame(Switzerland5, FG_matrix)
  
## Create the additive species variables and add them to Switzerland5.
## Switzerland5 will contain the new variables: p1_add, p2_add, p3_add and p4_add.
  ADD_matrix <- DI_data(prop = c("p1","p2","p3","p4"), data = Switzerland, what = "ADD")
  Switzerland5 <- data.frame(Switzerland5, ADD_matrix)
  
## Create all pairwise interaction variables and add them to Switzerland5.
## Switzerland5 will contain the new variables: p1.p2, p1.p3, p1.p4, p2.p3, p2.p4, p3.p4.
  FULL_matrix <- DI_data(prop = c("p1","p2","p3","p4"), data = Switzerland, what = "FULL")
  Switzerland5 <- data.frame(Switzerland5, FULL_matrix)
  
################################
  

################################ 
  
#### Short worked example (as illustrated the Details section)
  
## Create a dataframe
  p1 <- c(0.1, 0.25)
  p2 <- c(0.2, 0.25)
  p3 <- c(0.3, 0.25)
  p4 <- c(0.4, 0.25)
  minidataset1 <- data.frame(p1,p2,p3,p4)
  
## Check the rows sum to 1
  rowSums(minidataset1[1:4]) 
  
## Create the FG variables, assume two functional groups and theta = 1
  FG_matrix <- DI_data(prop = c("p1","p2","p3","p4"), FG = c("G","G","L","L"),
                       data = minidataset1, what = "FG")
  minidataset2 <- data.frame(minidataset1, FG_matrix)
  
## Create the FG variables, assume two functional groups and theta = 0.5
  FG_matrix <- DI_data(prop = c("p1","p2","p3","p4"), FG = c("G","G","L","L"), 
                       data = minidataset1, what = "FG", theta = 0.5)
  minidataset3 <- data.frame(minidataset1, FG_matrix)
## Create the ADD variables, assume theta = 0.5
  ADD_matrix <- DI_data(prop = c("p1","p2","p3","p4"), 
                        data = minidataset1, what = "ADD", theta = 0.5)
  minidataset3 <- data.frame(minidataset3, ADD_matrix)
## Add "_theta" to the new interaction variables to differentiate from when theta = 1
  names(minidataset3)[5:11] <- paste0(names(minidataset3)[5:11], "_theta")   

  
################################
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
