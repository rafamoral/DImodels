\name{functional_redundancy}
\alias{fr_communities}
\alias{fr_test}
\alias{fr_test_DI}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Test if a set of species are functionally redundant
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~

These functions can be used to assess if a set of species are functionally redundant or not, based on equivalence tests (Lakens, 2017) between community predictions. Users have two workflows 
\enumerate{
  \item \code{fr_test_DI} should be used for models fit using \code{\link{DI}} or \code{\link{autoDI}} functions and needs only a fitted model object along with names of species to test for redundancy.
  \item \code{fr_test} can be used for models fit using any other modelling package such as \code{\link{lm}()}, \code{\link{glm}()}, \code{nlme()} etc. or if users desire greater flexiblilty.
  In this case, users first need to run \code{fr_communities()} to get the set communities needed for comparisons, add any extra columns required for model predictions and pass the data and model object to \code{fr_test()}.
}

}
\usage{
# Workflow for model fit using DI() or autoDI()
fr_test_DI(object, redundant, 
           delta = delta_sd(0.5), alpha.level = 0.05)

# Workflow for model fit using any other modelling package.
fr_communities(prop, redundant)
fr_test(data, object, int_cols, 
        delta = delta_sd(0.5), alpha.level = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A model containing species proportions and optionally, their interactions as well as any additional treatment variable as predictors. For \code{fr_test_DI()} this should be a model fitted using \code{\link{DI}} or \code{\link{autoDI}}, but for \code{fr_test()} it can be fit using any popular modelling package such as \code{\link{lm}()}, \code{\link{glm}()}, \code{nlme()} etc.}
  \item{prop}{A vector of column names identifying the species proportions in the data used for fitting the model. For example, if the species proportions columns used as predictors are labelled p1 to p4, then \code{prop = c("p1","p2","p3","p4")}. Not needed for \code{fr_test_DI()} and will be inferred automatically from the model.}
  \item{redundant}{A vector containing names of species to be tested for functional redundancy. Should have at least two species and names should match those specified in \code{prop}.}
  \item{data}{The data.frame containing species proportions returned by \code{\link{fr_communities}()} with all additional columns needed for generating model predictions appended. For example, users may add columns corresponding to species interactions or to assess functional redundancy at specific values of a treatment variable.}
  \item{int_cols}{A vector containing the names or positions of columns in \code{data} corresponding to species interactions added by the user. Required only when using \code{fr_test}. If there are no interactions in the model, specify as NULL.}
  \item{delta}{Equivalence margin to use for comparisons. Can be specified in one of the following forms: 
    \enumerate{
    \item A function call of the form \code{\link{delta_sd}(x)}, where \code{x} is a positive multiplier. The resulting equivalence margin is \code{+/- x} times the residual standard error of \code{object}.
    \item A single numeric value, in which case the equivalence margin is \code{(-value, value)}.
    \item A numeric vector of length two, \code{c(x1, x2)}, specifying the lower and upper margin boundaries, respectively. This can be used to specify an asymmetric equivalence margin.
  }
  Defaults to \code{\link{delta_sd}(0.5)}. See \emph{Choosing equivalence margin} section in \strong{Details} section for more information.}
  \item{alpha.level}{A value between 0 and 1, specifying the singificance level for the comparisons. Defaults to 0.05.}
}
\details{

\emph{Functional redundancy}

When two or more species perform similar roles in an ecosystem such that substituting one or more species with others has little to no impact on a particular ecosystem response, the species are said to be functionally redundant with respect to that ecosystem function (Loreau et al., 2000, Kirwan et al., 2009).

Kirwan et al., 2009 proposed a test for assessing functional redundancy between two or more species in an ecosystem using the Diversity-Interactions (DI) framework. For a set (\eqn{F}) of \eqn{k} species to be functionally redundant, they should fulfil the following three conditions

\enumerate{
\item{Have equal identity effects i.e, monoculture performances (\eqn{\beta_1 = \beta_2 = ... = \beta_k})}.
\item{Have no within-group interactions among themselves (i.e., \eqn{\delta_ij = 0; \forall{i,j \in F})}}).
\item{Have identical interaction magnitudes with all remaining species in the ecosystem (i.e., \eqn{\delta_{iq} = \delta_{jq}; \forall{i,j \in F \, \& \, q \in F'}}).}
}

These conditions can be tested by fitting a series of heirarchical models and comparing them using an omnibus F-test. 
However, in this setting functional redundancy is implicitly assumed as the null hypothesis, effectively testing for the absence of functional redundancy rather than its presence. This is problematic because failing to reject the null hypothesis cannot be interpreted as evidence of functional redundancy and may lead to false conclusions.
Furthermore, the tests could reliably detect functional redundancy only when \eqn{\theta = 1} and returned false positives when \eqn{\theta \ne 1}.

The \code{fr_test} and \code{fr_test_DI} functions solve these problems by using equivalence testing (Lakens, 2017) to formally test for the presence of functional redundancy among species. Given a model object and a set \eqn{F} containing \eqn{k} species to be tested for redundancy, the function performs a series of equivalence tests comparing the model predicted responses for the following species communities:

\enumerate{
\item{Are the monoculture performances of all \eqn{k} species in \eqn{F} equivalent?}
\item{Is the performance of the equi-proportional mixture containing all \eqn{k} species in \eqn{F} equivalent to their monoculture average (i.e., is overyielding equivalent to zero)?}
\item{For each species in \eqn{F}, are the performances of its equi-proportional binary mixtures with every other species in \eqn{F'} equivalent?}
}

These tests are conceptually identical to those proposed in Kirwan et al., 2009, but are formulated as equivalence tests, allowing the null hypothesis to correspond to the absence of functional redundancy.
Furthermore, working directly with predictions from the final selected model eliminates the need for a separate model selection procedure for functional redundancy as well as any oddities due to \eqn{\theta}.

As a short example, if an ecosystem had 6 species (called \eqn{Sp_1}, \eqn{Sp_2}, \eqn{Sp_3}, \eqn{Sp_4}, \eqn{Sp_5}, and \eqn{Sp_6}) and we wished to test if species \eqn{Sp_1}, \eqn{Sp_2}, and \eqn{Sp_3} are functionally redundant. The following equivalence tests would be conducted:
\enumerate{
\item{Do the monocultures of \eqn{Sp_1}, \eqn{Sp_2}, and \eqn{Sp_3} have equivalent performance?}
\item{Is the performance of the equi-proportional mixture containing \eqn{Sp_1}, \eqn{Sp_2}, and \eqn{Sp_3} equivalent to their monoculture average?}
\item{Are the performances of binary mixtures of \eqn{Sp_1} and \eqn{Sp_4}, \eqn{Sp_1} and \eqn{Sp_5}, \eqn{Sp_1} and \eqn{Sp_6} equivalent? Similarly for \eqn{Sp_2} and \eqn{Sp_3}?}
}

The following figure shows a conceptual diagram of the communities that will be compared to assess functional redundancy between \eqn{Sp_1}, \eqn{Sp_2}, and \eqn{Sp_3}. 
Pie-glyphs identity the specific communities with \eqn{Sp_1}, \eqn{Sp_2}, and \eqn{Sp_3} shown using shades of green while the remaining species shown in shades of purple.

\figure{FR_example.png}{
  Conceptual diagram showing the specific communities that will be compared to assess functional redundancy between species in an ecosystem.
}

\emph{Choosing equivalence margin (i.e., delta)}

There is no objective theoretical boundary for defining the equivalence margin. 
It is commonly recommended to choose an equivalence margin based on domain expertise
or the smallest effect size of interest (SESOI). If a SESOI is not readily available,
it can be informed by study constraints such as the sample size, desired significance 
(\eqn{\alpha}) and power (\eqn{1 - \beta}) levels. 

Lakens (2019) recommends to define the equivalence margin based on the smallest
effect size which can be reliably detected in a study. However, this is implicitly
dependent on the study sample size and may result in a very large equivalence margin 
for studies with smaller sample sizes. 

In the absence of a defined domain-specific equivalence margin,
a reasonable starting point is to define the margin using a standardised difference scale.
For example, an equivalence margin corresponding to a small, medium, or large effect size, 
quantified using Cohen's \eqn{d} of 0.2, 0.5, or 0.8, respectively, can be used. 
The default maring used in the function is set as \eqn{0.5 \times \sigma} where sigma is 
residual standard deviation of the model, corresponding to a moderate effect size on the 
standardised difference scale. However, this should be used as a last resort and users
should define a margin justified by domain expertise or study-specific considerations.

Our recommendations, in decreasing order of preference are:
\enumerate{
\item{Define delta using domain expertise. A good starting point may be to define it as a percentage of the commonly observed range for the response being tested. However, this should not simply be based on the observed range of data being tested but inferred from similar previous studies involving the response of interest.}
\item{Define delta based on the smallest effect size that can be detected in the study. Given the study sample size, desired significance (\eqn{\alpha}) and power (\eqn{1 - \beta}) levels, conduct a power analysis determine the smallest effect that the study can reliably detect and use this as a basis for the equivalence margin.}
\item{Define delta using a standardised effect-size benchmark. In the absence of the previous two options, conventional benchmarks corresponding to a small, medium, or large effect sizes may provide a reasonable starting point for an equivalence margin.}
}

Regardless of how the equivalence margin is determined, it should always be specified 
before conducting the test. Furthermore, the chosen margin should always be reported
along with the test results and p-values. This would allow readers to interpret the 
result in context of the chosen equivalence margin and assess if the observed result
is appropriate for the specific response. It is also recommended to perform
a sensitivity analysis across a range of delta values to assess how conclusions
depend on the chosen equivalence margin.

%%  ~~ If necessary, more details than the description above ~~
}
\value{

All functions return a \code{data.frame} but with different columns.

\code{fr_communities} returns a \code{data.frame} containing species proportions for the communities required to assess functional redundancy among species specified in \code{redundant}. 
\describe{
  \item{Set}{Has values "Monos", "Within", or "Between" corresponding to whether the community is a monoculture, equi-proportional mixture for species specified in \code{redundant}, or binary mixtures with species not in `redundant` respectively. Internally, pairwise contrasts are computed for each community within a Set.}
  \item{Comm_ID}{Unique ID for each community in data. Used internally.}
}

\code{fr_test} and \code{fr_test_DI} return a \code{data.frame} containing results of pairwise comparisons between communities with the following columns.

\describe{
  \item{Set}{Has values "Monos", "Within", or "Between" corresponding to monoculture, within group and between group comparisons.}
  \item{Equiv_group}{Similar to \code{Set} but for between group comparisons it identifies the particular species being compared.}
  \item{Comm1}{An identifier for community 1 in the contrast.}
  \item{Comm2}{An identifier for community 2 in the contrast.}
  \item{Estimate}{The contrast estimate.}
  \item{SE}{Standard error for the contrast.}
  \item{df}{Residual df from the model.}
  \item{alpha}{Significance level used for equivalence tests.}
  \item{CI_level}{Confidence interval (derived using 1-2*alpha.level) for equivalent tests.}
  \item{CI_lo}{Lower limit of confidence interval around \code{Estimate}.}
  \item{CI_hi}{Upper limit of confidence interval around \code{Estimate}.}
  \item{delta_lo}{Lower bound of equivalence margin.}
  \item{delta_hi}{Upper bound of equivalence margin.}
  \item{t_raw}{Regular t-statistic for the contrast, calculated as the \code{Estimate} divided by \code{SE}.}
  \item{p_raw}{p-value from the standard two-sided t-test comparing if contrast estimate significantly differs from zero at the \code{alpha} significance level.}
  \item{p_lower}{p-value for the one-sided t-test assessing whether \code{Estimate} is greater than \code{delta_lo}.}
  \item{p_upper}{p-value for the one-sided t-test assessing whether \code{Estimate} is smaller than \code{delta_hi}.}
  \item{p_tost}{Overall p-value for the equivalence test. Maximum between \code{p_lower} and \code{p_upper}.}
  \item{Equivalent}{A logical value indicating whether the pairwise contrast is equivalent to zero (TRUE) or not (FALSE).}
}

}
\references{
%% ~put references to the literature/web site here ~

Lakens, D. (2017). Equivalence tests: A practical primer for t-tests, correlations, and meta-analyses. Social Psychological and Personality Science, 8(4), 355-362. doi:10.1177/1948550617697177.

Loreau, M., 2000. Biodiversity and ecosystem functioning: recent theoretical advances. Oikos,
91(1), pp.3-17

Kirwan L, J Connolly, JA Finn, C Brophy, A \enc{Lüscher}{}, D Nyfeler and MT Sebastia (2009) Diversity-interaction modelling - estimating contributions of species identities and interactions to ecosystem function. Ecology, 90, 2032-2038. 

Harms, C. and Lakens, D., 2018. Making 'null effects' informative: statistical techniques and inferential frameworks. Journal of clinical and translational research, 3(Suppl 2), p.382.

Lakens, D. and Delacre, M., 2020. Equivalence testing and the second generation p-value. Meta-Psychology, 4.

}

\author{
Rishabh Vishwakarma, Matthias Suter, James O’Malley, John A. Finn, Rafael A. Moral, and Caroline Brophy
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{DI}}
\code{\link{autoDI}}
\code{\link{fr_helpers}}
\code{\link{compare_communities}}

%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

## Load the Switzerland data
  data(Switzerland)
## Summarise the Switzerland data
  summary(Switzerland)

## Fit a DI model
  m1 <- DI(y = "yield", prop = 4:7, treat = 'nitrogen', DImodel = 'AV',
  density = 'density', estimate_theta = TRUE, data = Switzerland)
  summary(m1)
  
## Test functional redundancy between species p1 and p2 
## use wrapper function for DI models
  fr_test_DI(object = m1, redundant = c("p1", "p2"))

## By default delta is set to be 0.5*residual SE of model, 
## but user is urged to choose their own delta and it should be
## chosen before performing the test

## A good starting choice could be to set delta as some percentage of
## range of commonly observed values for the response across several experiments. 
## For example, one could assume that DM yields range between 5 and 20 t/ha and 
## thus they would set delta as 10% of 20, i.e., they consider a difference of 
## 2 t/ha to be neglible.
  my_delta <- 20*0.1
  print(my_delta)
  t1 <- fr_test_DI(object = m1, redundant = c("p1", "p2"),
                   delta = my_delta)

## The overall test result along with results for failing criteria are printed
  print(t1)
## We fail to establish functional redundancy between p1 and p2 because their
## monoculture performances and p1:p2 interaction were not found to be equivalent

## The summary() function can be used to get p-values for each individual contrast
  summary(t1)

## A plot() function is also provided to visualise the test results
  plot(t1)


## Users can also run the function with models fit using other modelling packages
## Here's an example using a model fit with glm, but same workflow applies for 
## models fit with other packages
  
  # Custom interaction terms for the model (just an example, can be anything)
  new_data <- Switzerland
  new_data$int1 <- new_data$p1*new_data$p2
  new_data$int2 <- new_data$p3*(new_data$p1 + new_data$p2)
  new_data$int3 <- new_data$p4*(new_data$p1 + new_data$p2)
  
  # Model
  m2 <- glm(yield ~ 0 + p1 + p2 + p3 + p4 + int1 + int2 + int3 + nitrogen,
            data = new_data)
  summary(m2)
  
  # First create specific communties that need to be compared using fr_communities()
  comms <- fr_communities(
              # prop is names of species in the model
              prop = c("p1", "p2", "p3", "p4"), 
              # This time we'll test for redundancy between p3 and p4.
              redundant = c("p3", "p4"))
  
  # Next we add any missing predictors from the model to this data
  ## Interactions (with same structure as specified when fitting the model)
  comms$int1 <- comms$p1*comms$p2
  comms$int2 <- comms$p3*(comms$p1 + comms$p2)
  comms$int3 <- comms$p4*(comms$p1 + comms$p2)
  
  ## Nitrogen treatment (compare at 150N level)
  comms$nitrogen <- factor("150", levels = c("50", "150"))
  
  # Now we can conduct test for functional redundancy
  t2 <- fr_test(
    # Model object
    object = m2,
    # Data with species interactions and extra predictors added in 
    data = comms,
    # Name of columns in data containing species interactions
    int_cols = c("int1", "int2", "int3"),
    # As an example, let's set delta as -1 and -2 (asymmetric interval)
    delta = c(-1, 2),
    # Perform test at 1% significance level
    alpha.level = 0.01
  )
  
  # Viewing test results
  print(t2)
  summary(t2)
  plot(t2)
  
  # If users with to explore the underlying contrast matrix used 
  # for computations, it can be accessed using
  attr(t2, "contrast_matrix")


}
