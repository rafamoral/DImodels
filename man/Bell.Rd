\name{Bell}
\alias{Bell}
\encoding{UTF-8}
\docType{data}
\title{
The "Bell" Dataset
} 
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
This dataset comes from a bacterial biodiversity experiment (Bell et al 2005). The bacterial ecosystems used were from semi-permanent rainpools that form in bark-lined depressions near the base of large European beech trees (\emph{Fagus sylvatica}). Microcosms consisting of sterile beech leaf disks and 10 ml of liquid (phosphate buffer) were inoculated with random combinations of 72 bacterial species isolated from these ecosystems. A total of 1,374 microcosms were constructed at richness levels of 1, 2, 3, 4, 6, 8, 9, 12, 18, 24, 36 and 72 species. The daily respiration rate of the bacterial
community in each microcosm was measured over three time intervals (days 0-7, 7-14 and
14-28) and the average over the three time intervals was recorded.
}

\usage{data("Bell")}
\format{
  A data frame with 1374 observations on the following 76 variables:
  \describe{
    \item{\code{id}}{A numeric vector uniquely identifying each row of the dataset.}
    \item{\code{community}}{A numeric vector identifying each unique community, i.e., two rows with the same community value also share the same set of p1 to p72 values.}
    \item{\code{richness}}{The number of species included in the initial composition, i.e., the number of proportions from p1 to p72 that are >0.}
    \item{\code{p1}}{A numeric vector indicating the initial proportion of species 1 in the community.}
    \item{\code{p2}}{A numeric vector indicating the initial proportion of species 2 in the community.}
    \item{\code{p3}}{A numeric vector indicating the initial proportion of species 3 in the community.}
    \item{\code{p4}}{A numeric vector indicating the initial proportion of species 4 in the community.}
    \item{\code{p5}}{A numeric vector indicating the initial proportion of species 5 in the community.}
    \item{\code{p6}}{A numeric vector indicating the initial proportion of species 6 in the community.}
    \item{\code{p7}}{A numeric vector indicating the initial proportion of species 7 in the community.}
    \item{\code{p8}}{A numeric vector indicating the initial proportion of species 8 in the community.}
    \item{\code{p9}}{A numeric vector indicating the initial proportion of species 9 in the community.}
    \item{\code{p10}}{A numeric vector indicating the initial proportion of species 10 in the community.}
    \item{\code{p11}}{A numeric vector indicating the initial proportion of species 11 in the community.}
    \item{\code{p12}}{A numeric vector indicating the initial proportion of species 12 in the community.}
    \item{\code{p13}}{A numeric vector indicating the initial proportion of species 13 in the community.}
    \item{\code{p14}}{A numeric vector indicating the initial proportion of species 14 in the community.}
    \item{\code{p15}}{A numeric vector indicating the initial proportion of species 15 in the community.}
    \item{\code{p16}}{A numeric vector indicating the initial proportion of species 16 in the community.}
    \item{\code{p17}}{A numeric vector indicating the initial proportion of species 17 in the community.}
    \item{\code{p18}}{A numeric vector indicating the initial proportion of species 18 in the community.}
    \item{\code{p19}}{A numeric vector indicating the initial proportion of species 19 in the community.}
    \item{\code{p20}}{A numeric vector indicating the initial proportion of species 20 in the community.}
    \item{\code{p21}}{A numeric vector indicating the initial proportion of species 21 in the community.}
    \item{\code{p22}}{A numeric vector indicating the initial proportion of species 22 in the community.}
    \item{\code{p23}}{A numeric vector indicating the initial proportion of species 23 in the community.}
    \item{\code{p24}}{A numeric vector indicating the initial proportion of species 24 in the community.}
    \item{\code{p25}}{A numeric vector indicating the initial proportion of species 25 in the community.}
    \item{\code{p26}}{A numeric vector indicating the initial proportion of species 26 in the community.}
    \item{\code{p27}}{A numeric vector indicating the initial proportion of species 27 in the community.}
    \item{\code{p28}}{A numeric vector indicating the initial proportion of species 28 in the community.}
    \item{\code{p29}}{A numeric vector indicating the initial proportion of species 29 in the community.}
    \item{\code{p30}}{A numeric vector indicating the initial proportion of species 30 in the community.}
    \item{\code{p31}}{A numeric vector indicating the initial proportion of species 31 in the community.}
    \item{\code{p32}}{A numeric vector indicating the initial proportion of species 32 in the community.}
    \item{\code{p33}}{A numeric vector indicating the initial proportion of species 33 in the community.}
    \item{\code{p34}}{A numeric vector indicating the initial proportion of species 34 in the community.}
    \item{\code{p35}}{A numeric vector indicating the initial proportion of species 35 in the community.}
    \item{\code{p36}}{A numeric vector indicating the initial proportion of species 36 in the community.}
    \item{\code{p37}}{A numeric vector indicating the initial proportion of species 37 in the community.}
    \item{\code{p38}}{A numeric vector indicating the initial proportion of species 38 in the community.}
    \item{\code{p39}}{A numeric vector indicating the initial proportion of species 39 in the community.}
    \item{\code{p40}}{A numeric vector indicating the initial proportion of species 40 in the community.}
    \item{\code{p41}}{A numeric vector indicating the initial proportion of species 41 in the community.}
    \item{\code{p42}}{A numeric vector indicating the initial proportion of species 42 in the community.}
    \item{\code{p43}}{A numeric vector indicating the initial proportion of species 43 in the community.}
    \item{\code{p44}}{A numeric vector indicating the initial proportion of species 44 in the community.}
    \item{\code{p45}}{A numeric vector indicating the initial proportion of species 45 in the community.}
    \item{\code{p46}}{A numeric vector indicating the initial proportion of species 46 in the community.}
    \item{\code{p47}}{A numeric vector indicating the initial proportion of species 47 in the community.}
    \item{\code{p48}}{A numeric vector indicating the initial proportion of species 48 in the community.}
    \item{\code{p49}}{A numeric vector indicating the initial proportion of species 49 in the community.}
    \item{\code{p50}}{A numeric vector indicating the initial proportion of species 50 in the community.}
    \item{\code{p51}}{A numeric vector indicating the initial proportion of species 51 in the community.}
    \item{\code{p52}}{A numeric vector indicating the initial proportion of species 52 in the community.}
    \item{\code{p53}}{A numeric vector indicating the initial proportion of species 53 in the community.}
    \item{\code{p54}}{A numeric vector indicating the initial proportion of species 54 in the community.}
    \item{\code{p55}}{A numeric vector indicating the initial proportion of species 55 in the community.}
    \item{\code{p56}}{A numeric vector indicating the initial proportion of species 56 in the community.}
    \item{\code{p57}}{A numeric vector indicating the initial proportion of species 57 in the community.}
    \item{\code{p58}}{A numeric vector indicating the initial proportion of species 58 in the community.}
    \item{\code{p59}}{A numeric vector indicating the initial proportion of species 59 in the community.}
    \item{\code{p60}}{A numeric vector indicating the initial proportion of species 60 in the community.}
    \item{\code{p61}}{A numeric vector indicating the initial proportion of species 61 in the community.}
    \item{\code{p62}}{A numeric vector indicating the initial proportion of species 62 in the community.}
    \item{\code{p63}}{A numeric vector indicating the initial proportion of species 63 in the community.}
    \item{\code{p64}}{A numeric vector indicating the initial proportion of species 64 in the community.}
    \item{\code{p65}}{A numeric vector indicating the initial proportion of species 65 in the community.}
    \item{\code{p66}}{A numeric vector indicating the initial proportion of species 66 in the community.}
    \item{\code{p67}}{A numeric vector indicating the initial proportion of species 67 in the community.}
    \item{\code{p68}}{A numeric vector indicating the initial proportion of species 68 in the community.}
    \item{\code{p69}}{A numeric vector indicating the initial proportion of species 69 in the community.}
    \item{\code{p70}}{A numeric vector indicating the initial proportion of species 70 in the community.}
    \item{\code{p71}}{A numeric vector indicating the initial proportion of species 71 in the community.}
    \item{\code{p72}}{A numeric vector indicating the initial proportion of species 72 in the community.}
    \item{\code{response}}{A numeric vector giving the average daily respiration rate of the bacterial community.} 
  }
}


\details{
%%  ~~ If necessary, more details than the __description__ above ~~
 
\strong{What are Diversity-Interactions (DI) models?} 
 

Diversity-Interactions (DI) models (Kirwan et al 2009) are a set of tools for analysing and interpreting data from experiments that explore the effects of species diversity on community-level responses. We strongly recommend that users read the short introduction to Diversity-Interactions models (available at: \code{\link{DImodels}}). Further information on Diversity-Interactions models is also available in Kirwan et al 2009 and Connolly et al 2013.

The \code{Bell} dataset is analysed using Diversity-Interactions models in both Brophy et al 2017 and Connolly et al 2013. 
 

}

\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
Bell T, JA Newman, BW Silverman, SL Turner and AK Lilley (2005) The contribution of species richness and composition to bacterial services. Nature, 436, 1157-1160.  
 
}
\references{
%references
Brophy C, A Dooley, L Kirwan, JA Finn, J McDonnell, T Bell, MW Cadotte and J Connolly (2017) Biodiversity and ecosystem function: Making sense of numerous species interactions in multi-species communities. Ecology, 98, 1771-1778. 
 
Connolly J, T Bell, T Bolger, C Brophy, T Carnus, JA Finn, L Kirwan, F Isbell, J Levine, A \enc{Lüscher}{}, V Picasso, C Roscher, MT Sebastia, M Suter and A Weigelt (2013) An improved model to predict the effects of changing biodiversity levels on ecosystem function. Journal of Ecology, 101, 344-355. 
 
Kirwan L, J Connolly, JA Finn, C Brophy, A \enc{Lüscher}{}, D Nyfeler and MT Sebastia (2009) Diversity-interaction modelling - estimating contributions of species identities and interactions to ecosystem function. Ecology, 90, 2032-2038. 
}
\examples{

## Load the Bell data
  data(Bell)
## View the first five entries
  head(Bell)
## Explore the variabes in sim1
  str(Bell)

## Check that the proportions sum to 1 (required for DI models)
## p1 to p72 are in the 4th to 75th columns in Bell
  Bellsums <- rowSums(Bell[4:75])
  summary(Bellsums)

## Check characteristics of Bell
  hist(Bell$response)
  summary(Bell$response)
  plot(Bell$richness, Bell$response)


\donttest{
## This code takes around 11 seconds to run

## Fit the average pairwise model using DI and the AV tag, with theta estimated
  m1 <- DI(y = "response", prop = 4:75, DImodel = "AV", estimate_theta = TRUE, data = Bell)
  summary(m1)
  CI_95 <- theta_CI(m1, conf = .95)
  CI_95
  plot(m1)
  library(hnp)
  hnp(m1)
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

}
\keyword{datasets}
