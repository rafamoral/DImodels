## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval = FALSE------------------------------------------------------------
#  install.packages("DImodels")
#  library("DImodels")

## ---- eval = FALSE------------------------------------------------------------
#  ?DImodels

## ---- eval = FALSE------------------------------------------------------------
#  ?sim3

## ---- echo = FALSE, results='asis'--------------------------------------------
library(DImodels)
data("design_a")
knitr::kable(head(design_a))

## ---- echo = TRUE, results='asis'---------------------------------------------
data("sim3")
knitr::kable(head(sim3, 10))

## ---- echo = TRUE-------------------------------------------------------------
hist(sim3$response, xlab = "Response", main = "")
# Similar graphs can also be generated for the other species proportions.
plot(sim3$p1, sim3$response, xlab = "Proportion of species 1", ylab = "Response")
summary(sim3$response)

## ---- echo = TRUE-------------------------------------------------------------
auto1 <- autoDI(y = "response", prop = 4:12, treat = "treatment", 
                FG = c("FG1","FG1","FG1","FG1","FG1","FG2","FG2","FG3","FG3"), data = sim3, 
                selection = "Ftest")

## ---- eval = FALSE------------------------------------------------------------
#  ?autoDI

## ---- echo = TRUE-------------------------------------------------------------
summary(auto1)

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  theta_CI(auto1, conf = .95)

## ---- echo = TRUE-------------------------------------------------------------
m1 <- DI(y = "response", prop = 4:12, 
         FG = c("FG1","FG1","FG1","FG1","FG1","FG2","FG2","FG3","FG3"), treat = "treatment", 
         DImodel = "FG", data = sim3)
summary(m1)

## ---- echo = TRUE-------------------------------------------------------------
m1_theta <- update_DI(object = m1, estimate_theta = TRUE)
coef(m1_theta)

## ---- echo = TRUE-------------------------------------------------------------
m2 <- DI(y = "response", prop = 4:12, 
         FG = c("FG1","FG1","FG1","FG1","FG1","FG2","FG2","FG3","FG3"), treat = "treatment", 
         DImodel = "FG", extra_formula = ~ (p1 + p2 + p3 + p4):treatment,
         data = sim3)
summary(m2)

## ---- echo = TRUE-------------------------------------------------------------
FG_matrix <- DI_data(prop = 4:12, FG = c("FG1","FG1","FG1","FG1","FG1","FG2","FG2","FG3","FG3"), 
                     data = sim3, what = "FG")
sim3a <- data.frame(sim3, FG_matrix)

## ---- echo = TRUE-------------------------------------------------------------
m3 <- DI(y = "response", prop = 4:12, 
         FG = c("FG1","FG1","FG1","FG1","FG1","FG2","FG2","FG3","FG3"),
         treat = "treatment", DImodel = "FG", 
         extra_formula = ~ (bfg_FG1_FG2 + bfg_FG1_FG3 + bfg_FG2_FG3 +
                              wfg_FG1 + wfg_FG2 + wfg_FG3) : treatment, data = sim3a)
summary(m3)

## ---- echo = TRUE-------------------------------------------------------------
sim3a$treatmentA <- as.numeric(sim3a$treatment == "A")

## ---- echo = TRUE-------------------------------------------------------------
m3 <- DI(y = "response",
         custom_formula = response ~ 0 + p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 +
           treatmentA + bfg_FG1_FG2 + bfg_FG1_FG3 + bfg_FG2_FG3, data = sim3a)
summary(m3)

