richness_vs_DI <- function(y, prop, data, extra_formula) {
  ## general set up
  if(missing(y)) {
    stop("You must supply a response variable name or column index through the argument 'y'.\n")
  }
  
  n <- nrow(data)
  if(missing(extra_formula)) {
    extra_formula <- 0
  }
  
  ## setting up dataset
  newdata <- data
  ## creating AV variable with theta = 0.5
  newdata$AV0.5 <- DI_data(prop = prop, data = data, what = "AV", theta = 0.5)
  ## creating richness variable
  newdata$richness <- apply(newdata[,prop], 1, function(x) sum(x > 0))
  ## defining response variable
  if(!is.character(y)) y <- names(data)[y]
  
  ## setting up extra formula (if present)
  if(!inherits(extra_formula, "formula")) {
    extra_terms <- ""
  } else {
    extra_terms <- paste("+", paste(extra_formula)[2])
  }
  
  ## fit the models
  ## model 1: the richness model
  formula_richness <- formula(paste(y, "~ richness", extra_terms))
  m1 <- lm(formula_richness, data = newdata)
  
  ## model 2: the DI equivalent
  formula_DI1 <- formula(paste(y, "~ AV0.5", extra_terms))
  m2 <- lm(formula_DI1, data = newdata)
  
  ## model 3: allowing different identity effects
  m3 <- suppressWarnings(
    suppressMessages(
      DI(y = y, prop = prop, DImodel = "AV", data = data, theta = 0.5,
         extra_formula = extra_formula)
      ))
  m3$df.residual <- m3$df.residual + 1
  
  ## model 4: estimating theta
  m4 <- suppressWarnings(
    suppressMessages(
      DI(y = y, prop = prop, DImodel = "AV", data = data, estimate_theta = TRUE,
         extra_formula = extra_formula)
    ))
  
  ## get AIC values for the four models
  all_AICs <- c(AIC(m1), AIC(m2), AIC(m3), AIC(m4))
  model_list <- list(m1, m2, m3, m4)
  
  ## print AIC values
  the_table <- data.frame("AIC" = all_AICs,
                          "df" = n - c(m1$df.residual,
                                       m2$df.residual,
                                       m3$df.residual,
                                       m4$df.residual),
                          "Description" = c("Model 1: Richness only",
                                            "Model 2: Average interactions 'AV' DImodel with common identity effects and theta = 0.5",
                                            "Model 3: Average interactions 'AV' DImodel with unique identity effects and theta = 0.5",
                                            "Model 4: Average interactions 'AV' DImodel with unique identity effects and theta estimated"),
                          row.names = NULL)
  
  message("\n", strrep("-", getOption("width")))
  message("\nInvestigating richness model and three DI alternatives\n")
  print(the_table, right = FALSE)
  
  ## return the model with smallest AIC
  message("\n", strrep("-", getOption("width")))
  message("richness_vs_DI is limited in terms of model selection. Only four models are explored. See ?autoDI and ?DI for more options.")
  message(strrep("-", getOption("width")))
  return(invisible(model_list[[which.min(all_AICs)]]))
}