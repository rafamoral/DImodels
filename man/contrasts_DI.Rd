\name{contrasts_DI}
\alias{contrasts_DI}
\alias{contrast_matrix}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Compute contrasts for Diversity Interactions (DI) models
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
\code{contrast_vars} calculates and tests contrasts using \code{\link[multcomp]{glht}} for a model object created by the \code{\link{DI}} or \code{\link{autoDI}} functions.

\code{contrast_matrix} is a helper function which can be used to prepare a contrast matrix which contains all the ID effects, interactions and additional model variables. The output of this function can then be passed to \code{contrasts_DI} to test contrasts.
}
\usage{
contrasts_DI(object, contrast_vars, contrast, verbose, ...)

contrast_matrix(object, contrast_vars)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A model output which is the output of \code{\link{DI}} or \code{\link{autoDI}}.}
  \item{contrast_vars}{A matrix or data-frame containing the values of the species prportions in a contrast style (i.e., positive and negative values). The species ID effects and interactions would be calculated automatically (the value of the non-linear parameter theta would also be incorporated) and added to the contrast. 
  
  The user can also specify any additional variable (like treatment, density) from the model to calculate the contrast at the specific value of that variable(s).
  
  Furthermore, it is also possible to manually specify the ID or interaction effects to have a specific value and override those calculated by the function.
  
  Any species proportion or additional variable not specified would be assumed to be 0. See examples showing the uses of these options.

  Will be overridden if \code{contrast} is specified.}
  \item{contrast}{A data-frame or matrix containing the values for generating the contrast. The number of columns in the matrix/data-frame should be the same as the number of linear coefficients in the model. It is assumed that the columns in \code{contrast} are in the same order as the coefficeints in the model and that the user would calculate the appropriate species interaction values by incorporating the theta value (if it's present) from the model. The matrix/data-frame can also be given row names to track the different contrasts. 
  
                  Recommended only for advanced users who want more control for generating the contrast as most general contrasts can be generated using the \code{contrast_vars} argument. 
    
                  Overrides \code{contrast_vars} if specified.}
  \item{verbose}{A boolean (TRUE/FALSE) value indicating whether to print the internally generated contrast matrix. Default is TRUE}
  \item{...}{Additional arguments passed to the \code{\link[multcomp]{glht}} function}
}
\details{
The contrasts are calculated and tested using the \code{\link[multcomp]{glht}} function in the \code{multcomp} package. If the various contrasts are not explicitly named, they'll be given names `Test 1`, `Test 2`, ..., `Test n`.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An object of class \code{\link[multcomp]{glht}} is returned
}
\references{
%% ~put references to the literature/web site here ~
Kirwan L, J Connolly, JA Finn, C Brophy, A \enc{Lüscher}{}, D Nyfeler and MT Sebastia (2009) Diversity-interaction modelling - estimating contributions of species identities and interactions to ecosystem function. Ecology, 90, 2032-2038. 

Torsten Hothorn, Frank Bretz and Peter Westfall (2008).
  Simultaneous Inference in General Parametric Models.
  Biometrical Journal 50(3), 346--363.
}
\author{
Rafael A. Moral, John Connolly, Rishabh Vishwakarma and Caroline Brophy
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{DI}}
\code{\link{autoDI}}
\code{\link[multcomp]{glht}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

## Load the Switzerland data
  data(Switzerland)
## Summarise the Switzerland data
  summary(Switzerland)

## Fit a DI model
  m1 <- DI(y = "yield", prop = 4:7, treat = 'nitrogen', DImodel = 'AV',
  density = 'density', estimate_theta = FALSE, data = Switzerland)
  summary(m1)
  
  
## Create contrasts 
## We test the following five sample contrasts
## 1. p1 monoculture against p2 monoculture
## 2. p1 and p2 binary (50-50) mixture against the p4 monoculture
## 3. p1 and p2 binary (50-50) mixture against p3 and p4 (50-50) binary mixture
## 4. Three species equi-proportional mixture (1/3 each) of p2, p3, and p4 against p1 monoculture
## 5. p4 monoculture against a 70-30 mixture of p2 and p3 respectively.
  contrast_vars <- data.frame("p1" = c(1,  0.5,  0.5,  -1,    0),
                              "p2" = c(-1, 0.5,  0.5, 1/3, -0.7),
                              "p3" = c(0,    0, -0.5, 1/3, -0.3),
                              "p4" = c(0,   -1, -0.5, 1/3,    1))
  con1 <- contrasts_DI(object = m1,
                       contrast_vars = contrast_vars)
                       
## Calling summary on the returned object would show a significane table for the tests
  summary(con1)                          

## In the previous example, note how the species ID effects and interaction term 
## was added to the contrast matrix. Additionally, the nitrogen50 and densityhigh
## variables were added with a value 0

## We could also calculate the contrast at a specific value for these additional variables
## Calculate contrast at nitrogen = 50 for the p1 monoculture against p3 monoculture
  contrast_vars2 <- data.frame("p1" = 1, "p3" = -1, "nitrogen50" = 1)
  # Name contrast to track it
  rownames(contrast_vars2) <- "p1 mono vs p3 mono at 50 N"
  con2 <- contrasts_DI(object = m1,
                       contrast_vars = contrast_vars2)
## Notice that p2 and p4 were added to the contrast with a value 0 
  summary(con2)                          

## It is also possible to specify the ID or interactions effects to have a set value
  contrast_vars3 <- data.frame("p2_ID" = 0.5, 
                               "p3_ID" = -0.5, 
                               "AV" = 0.375)
  con3 <- contrasts_DI(object = m1,
                       contrast_vars = contrast_vars3)
## Notice the values for terms specified are preserved while others are calculated or set to 0
  summary(con3)                          

#############################################################################
## Using the `contrast` agrument for creating contrasts
## Contrasts for difference between monocultures of p1 and p2, p3 and p4, p2 and p3
  con4 <- contrasts_DI(object = m1,
                       contrast = data.frame('p1vp2 Mono' = c(1, -1, 0, 0, 0, 0, 0),
                                             'p3vp4 Mono' = c(0, 0, 1, -1, 0, 0, 0),
                                             'p2vp3 Mono' = c(0, -1, 1, 0, 0, 0, 0)))
  summary(con4)
  
## Contrasts for 50:50 mixture of p1 and p2 vs 50:50 mixture of p3 and p4  
  con5 <- contrasts_DI(object = m1,
                       contrast = data.frame('p1p2 vs p3p4' = c(0.5, 0.5, -0.5, -0.5, 0, 0, 0)))
  summary(con5)
  
## There is also a helper function called `contrast_matrix` to create the 
## full contrast matrix with interaction terms and any additional variables
## in the model which can be modified further by the user and passed
## into the contrast parameter

## Compare the three species mixture of p1, p2 and p3 to the centroid mixture
  mix1 <- contrast_matrix(object = m1, 
                          contrast_vars = data.frame("p1" = 1/3, "p2" = 1/3,
                                                     "p3" = 1/3))
  mix2 <- contrast_matrix(object = m1, 
                          contrast_vars = data.frame("p1" = 1/4, "p2" = 1/4,
                                                     "p3" = 1/4, "p4" = 1/4))

## The interaction terms and nitrogen and density terms are all added
## Subtract the two vectors to get a contrast
  contr <- mix1 - mix2
  con6 <- contrasts_DI(object = m1,
                       contrast = contr)
  summary(con6)
  
## Could also modify a variable the returned output to test another contrast
## Suppose we wish to compare the three species mixture at 150 kg N vs the 
## centroid at 50 kg N
  mix1[, "nitrogen50"] <- 0
  mix2[, "nitrogen50"] <- 1
  contr <- mix1 - mix2
  con7 <- contrasts_DI(object = m1,
                       contrast = contr)
  summary(con7)
}
