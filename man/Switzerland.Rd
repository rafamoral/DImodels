\name{Switzerland}
\alias{Switzerland}
\encoding{UTF-8}
\docType{data}
\title{
The "Switzerland" Dataset
} 
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
This dataset comes from a grassland biodiversity experiment that was conducted in Switzerland as part of the "Agrodiversity Experiment" (Kirwan et al 2014). A total of 68 grassland plots were established across a gradient of species diversity, and two additional treatments (nitrogen fertiliser and total seed density) were also manipulated. The proportions of four species were varied across the plots: there were plots with 100\% of a single species, and 2- and 4-species mixtures with  varying proportions (e.g., (0.5, 0.5, 0, 0) and (0.7, 0.1, 0.1, 0.1)). Nitrogen fertiliser was either 50 or 100 kg N per annum and total seed density was either low or high. Total annual yield per plot was recorded for the first year after establishment. An analysis of the Switzerland dataset is presented in Kirwan et al 2009. 
}
\usage{data("Switzerland")}
\format{
  A data frame with 68 observations on the following 8 variables:
  \describe{
    \item{\code{plot}}{A numeric vector uniquely identifying each of the 68 plots.}
    \item{\code{nitrogen}}{A factor with two levels: "50" or "150" to indicate the level of nitrogen fertiliser (kg N per annum) applied to the plot.}
    \item{\code{density}}{A factor with two levels: "low" and "high" to indicate the level of total seed density used when sowing the plot.}
    \item{\code{p1}}{A numeric vector indicating the proportion of species 1 in the plot. Species 1 was the grass species \emph{Lolium perenne}.}
    \item{\code{p2}}{A numeric vector indicating the proportion of species 2 in the plot. Species 2 was the grass species \emph{Dactylis glomerata}.}
    \item{\code{p3}}{A numeric vector indicating the proportion of species 3 in the plot. Species 3 was the legume species \emph{Trifolium pratense}.}
    \item{\code{p4}}{A numeric vector indicating the proportion of species 4 in the plot. Species 4 was the legume species \emph{Trifolium repens}.}
    \item{\code{yield}}{A numeric vector giving the total dry matter yield for the plot (tonnes per hectare per annum).} 
  }
}


\details{
%%  ~~ If necessary, more details than the __description__ above ~~
\strong{What are Diversity-Interactions (DI) models?}

Diversity-Interactions (DI) models (Kirwan et al 2009) are a set of tools for analysing and interpreting data from experiments that explore the effects of species diversity on community-level responses. We strongly recommend that users read the short introduction to Diversity-Interactions models (available at: \code{\link{DImodels}}). Further information on Diversity-Interactions models is also available in Kirwan et al 2009 and Connolly et al 2013.

\strong{Functional groups}

In Ecology, species can be categorised into 'functional groups' based on their traits and functions. Here, the four species comprise two grasses (species 1 and 2) and two legumes (species 3 and 4); this is one possible 'functional group' categorisation of the four species. 

}

\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
Kirwan L, J Connolly, C Brophy, O Baadshaug, G Belanger, A Black, T Carnus, R Collins, J \enc{Čop}{}, I Delgado, A De Vliegher, A Elgersma, B Frankow-Lindberg, P Golinski, P Grieu, AM Gustavsson, \enc{Á}{} \enc{Helgadóttir}{}, M \enc{Höglind}{}, O Huguenin-Elie, M \enc{Jørgensen}{}, \enc{Ž}{} \enc{Kadžiuliene}{}, T Lunnan, A \enc{Lüscher}{}, P Kurki, C Porqueddu, MT Sebastia, U Thumm, D Walmsley and JA Finn (2014) The Agrodiversity Experiment: three years of data from a multisite study in intensively managed grasslands. Ecology, 95, 2680.
% ė character in "Kadžiulienė" generating errors in LaTeX manual
}
\references{
Connolly J, T Bell, T Bolger, C Brophy, T Carnus, JA Finn, L Kirwan, F Isbell, J Levine, A \enc{Lüscher}{}, V Picasso, C Roscher, MT Sebastia, M Suter and A Weigelt (2013) An improved model to predict the effects of changing biodiversity levels on ecosystem function. Journal of Ecology, 101, 344-355.

Kirwan L, J Connolly, JA Finn, C Brophy, A \enc{Lüscher}{}, D Nyfeler and MT Sebastia (2009) Diversity-interaction modelling - estimating contributions of species identities and interactions to ecosystem function. Ecology, 90, 2032-2038.
}
\examples{

## Load the Switzerland data
  data(Switzerland)
## View the first few entries
  head(Switzerland)
## Explore the variables in Switzerland
  str(Switzerland)

## Histogram of the response variable yield
  hist(Switzerland$yield)
## Explore the marginal relationship between yield and each predictor
  plot(Switzerland$p1, Switzerland$yield)
  plot(Switzerland$p2, Switzerland$yield)
  plot(Switzerland$p3, Switzerland$yield)
  plot(Switzerland$p4, Switzerland$yield)
  boxplot(yield ~ nitrogen, data = Switzerland)
  boxplot(yield ~ density, data = Switzerland)
  
## Check that the proportions sum to 1 (required for DI models)
## p1 to p4 are in the 4th to 7th columns in Switzerland
  Switzerlandsums <- rowSums(Switzerland[4:7])
  summary(Switzerlandsums)
  
## Model selection by F-test
  auto1 <- autoDI(y = "yield", density = "density", prop = c("p1","p2","p3","p4"), 
                  treat = "nitrogen", FG = c("G","G","L","L"), data = Switzerland, 
                  selection = "Ftest")
  summary(auto1)
  
## Fit the model chosen by autoDI using DI
  m1 <- DI(y = "yield", density = "density", prop = 4:7, DImodel = "FG", FG = c("G","G","L","L"), 
           data = Switzerland)
  summary(m1)  
  plot(m1)
  
\donttest{
## Check goodness-of-fit using a half-normal plot with a simulated envelope
  library(hnp)
  hnp(m1)
}
  
## Set up the functional group interactions and add to a new Switzerland2 dataset
  FG_matrix <- DI_data(prop = 4:7, FG = c("G","G","L","L"), 
                          data = Switzerland, what = "FG")
  Switzerland2 <- data.frame(Switzerland, FG_matrix)
  
## Additional model testing using DI to test for interactions with nitrogen
  m2 <- DI(y = "yield", block = "density", prop = 4:7, DImodel = "FG", FG = c("G","G","L","L"), 
           data = Switzerland2, extra_formula = ~ nitrogen:bfg_G_L)
  summary(m2) 

}
\keyword{datasets}
