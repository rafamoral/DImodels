\name{compare_communities}
\alias{compare_communities}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Compare communities and assign compact letter displays (CLD)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
\code{compare_communities} performs comparisons between species and assigns each a compact letter display (CLD) using \code{\link[multcomp]{glht}} and \code{\link[multcompView]{multcompLetters}} for a model object created by the \code{\link{DI}} or \code{\link{autoDI}} functions.

}
\usage{
compare_communities(object, data, ref = NULL, 
                    adjust = formals(multcomp::adjusted)$type, 
                    verbose = FALSE, alpha.level = 0.05,
                    Letters = c(letters, LETTERS, "."),
                    ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A model output which is the output of \code{\link{DI}} or \code{\link{autoDI}}.}
  \item{data}{A data.frame containing the proportions for the species communities which are to be compared. Any additional variables such as specific values for a treatment variable can also be specified.}
  \item{ref}{An optional argument specifying a particular to be considered a reference. If specified, all comparisons would be made with respect to this community, as opposed to pairwise comparisons.}
  \item{adjust}{A character string specifying the p-value adjustment method to be used to account for multiple comparisons. All arguments supported by the \code{\link[multcomp]{adjusted}} are supported. The default is "single-step" adjustment.}
  \item{verbose}{A boolean (TRUE/FALSE) value indicating whether to print the internally generated contrast matrix. Default is TRUE}
  \item{alpha.level}{A value between 0 and 1, specifying the singificance level for the comparisons. Defaults to 0.05.}
  \item{Letters}{The characters to be assigned to the communities to indicate significant differences. Defaults to lowercase latin characters (i.e., a, b, c, etc.).}
  \item{...}{Additional arguments passed to the \code{\link[multcomp]{glht}}, \code{\link[multcomp]{adjusted}}, or \code{\link[multcompView]{multcompLetters}} functions}
}
\details{
The contrasts are calculated and tested using the \code{\link[multcomp]{glht}} function in the \code{multcomp} package while the cld letters are assigned using the \code{\link[multcompView]{multcompLetters}} function from \code{multcompView}.
The rownames of the specified \code{data.frame} are used to identify the communities and are labelled 1, 2, ..., n, etc. by default, but this can be changed by updating the rownames of the specified \code{data.frame}.

Note: If any variables needed for making predictions are missing in `data`, the would be assumed to be at their median value (if numeric) or base level (if categorical). Specify the values for all variables manually to perform comparisons at a specific level of any categorical variable (treatment for example).
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An list containing the following objects is returned.

\item{Contrasts}{An object of class \code{\link[multcomp]{summary.glht}} containing the summary of the specific comparisons performed internally for assigning the compact display letters.}
\item{CLD}{An object of class \code{\link[multcompView]{multcompLetters}} with compact letters assigned to each species community based on the p-values shown in the `Contrasts` object.}

}
\references{
%% ~put references to the literature/web site here ~

Piepho, Hans-Peter (2004) "An Algorithm for a Letter-Based Representation of All-Pairwise Comparisons", Journal of Computational and Graphical Statistics, 13(2)456-466.

Torsten Hothorn, Frank Bretz and Peter Westfall (2008).
  Simultaneous Inference in General Parametric Models.
  Biometrical Journal 50(3), 346--363.
  

Frank Bretz, Torsten Hothorn and Peter Westfall (2010), Multiple Comparisons Using R, CRC Press, Boca Raton.

Jason C. Hsu (1996), Multiple Comparisons. Chapman & Hall, London.

Kirwan L, J Connolly, JA Finn, C Brophy, A \enc{Lüscher}{}, D Nyfeler and MT Sebastia (2009) Diversity-interaction modelling - estimating contributions of species identities and interactions to ecosystem function. Ecology, 90, 2032-2038. 
}

\author{
Rafael A. Moral, John Connolly, Rishabh Vishwakarma and Caroline Brophy
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{DI}}
\code{\link{autoDI}}
\code{\link{contrasts_DI}}
\code{\link[multcomp]{glht}}
\code{\link[multcomp]{adjusted}}
\code{\link[multcompView]{multcompLetters}}

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
  
## Basic comparison of the first five observations in the Switzeraland data
  head(Switzerland, 5)
  compare_communities(m1, data = Switzerland[1:5, ])

## The pairiwise comparisons between each community as well as the cld letters are returned.

## By default, single-step p-value adjustment are performed to account for multiple comparisons. 
## This can be changed using adjust argument. For example "none" would perform no adjustment
  compare_communities(m1, data = Switzerland[1:5, ],
                      adjust = "none")

## Bonferroni adjustment
  compare_communities(m1, data = Switzerland[1:5, ],
                      adjust = "bonferroni")
                      
## If any variable needed for predictions is missing in the data, 
## a warning will be thrown and sensible values will be used.   
## In this example, we only specify the species proportions and nothing for nitrogen and density
  compare_communities(m1, data = Switzerland[1:5, 4:7],
                      adjust = "bonferroni")
                      
## The alpha.level and symbols to be used for cld can be changed as follows
  compare_communities(m1, data = Switzerland[1:5, ], 
                      adjust = "none", alpha.level = 0.01,
                      Letters = LETTERS)

## Change row.names of data to get informative labels for the tests
  comp_data <- Switzerland[1:5, ]
  rownames(comp_data) <- c("p1_dom", "p2_dom", "p3_dom", "p4_dom", "Centroid")
  
  compare_communities(m1, data = comp_data, 
                      adjust = "none", alpha.level = 0.01,
                      Letters = LETTERS)

## Comparisons can also be performed against a specific community as reference 
## using the `ref` argument. Using "Centroid" as ref in this example, 
## can also be specified using row-number, i.e., 5
  compare_communities(m1, data = comp_data, ref = 5,    
                      adjust = "none", 
                      Letters = LETTERS)

}
