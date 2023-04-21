# y = response variable index or name
# block = block variable index or name
# density = density variable index or name
# prop = vector of indices or names of the species proportions
# treat = environmental covariate index or name
# data = dataset in data.frame format
# selection = "Ftest" to perform F (or X2) tests to select best model,
#             "AIC" to use select model with smallest AIC

autoDI <- function(y, prop, data, block, density, treat, FG = NULL, 
                   selection = c("Ftest","AIC","AICc","BIC","BICc"), 
                   step0 = FALSE, step4 = TRUE) {
  if(missing(y)) stop("You must supply a response variable name or column index through the argument 'y'.\n")
  
  ########## RV change ###############
  if (inherits(data, 'tbl')){
    data <- as.data.frame(data)
  }
  
  ####################################
  
  # family / binomial denominator lock
  # include family and total as function arguments to lift the lock
  family <- "gaussian"
  total <- NULL
  # family lock
  if(!(family %in% c("gaussian","normal")))
    stop("As of version ", packageVersion("DImodels"),
         " DI models are implemented for family = 'gaussian' (= 'normal') only")
  
  # set theta flag
  estimate_theta <- TRUE
  # check if user chose "Ftest", "AIC", "AICc" or "BIC"
  selection <- match.arg(selection)
  # default family set as normal
  if(missing(family) || family == "normal") family <- "gaussian"
  # checks for binomial
  if(family %in% c("binomial","quasibinomial")) {
    if(missing(total)) {
      if(any(!(data[,y] %in% c(0,1)))) {
        stop("total must be informed for non-binary discrete proportion data")
      } else total <- 1
    } else total <- data[,total]
  } else total <- NULL
  # checks for quasi and information criteria
  if(family %in% c("quasipoisson","quasibinomial") & selection %in% c("AIC","AICc","BIC","BICc"))
    stop("cannot compute information criteria for quasi models, use selection = 'Ftest'")
  # flags if block/density and/or treat are missing
  if(missing(block)) block <- NA
  if(missing(density)) density <- NA
  treat_flag <- FALSE
  if(missing(treat)) {
    treat <- NA
    treat_flag <- TRUE
  }
  
  if(missing(FG)){
    FG <- NULL
  }
  
  # Step 0
  if(step0) {
    autoDI_step0(y = y, block = block, density = density, treat = treat, family = family, data = data)
  }
  
  # Step 1
  step1_model <- autoDI_step1(y = y, block = block, density = density, prop = prop, treat = treat, family = family, data = data, selection = selection)
  
  selected_model <- step1_model$model
  theta_flag <- step1_model$theta_flag
  
  # step 2
  step2_model <- autoDI_step2(y = y, block = block, density = density, prop = prop, treat = treat, family = family, FG = FG, data = data, theta_flag = theta_flag, selection = selection, selected_model = selected_model)
  
  # step 3
  step3_model <- autoDI_step3(selected_model = step2_model, selection = selection, family = family)
  
  # step 4
  # If lack of fit not possible skip it
  if(step4 & (nrow(unique(data[, prop])) == nrow(data))){
    message("\n", strrep("-", getOption("width")))
    message("Step 4: Comparing the final selected model with the reference (community) model")
    message('Lack of fit test is not possible as there are no repititions of communities in the data. Skipping step 4')
    step4 <- FALSE
  }
  
  if(step4) {
    autoDI_step4(prop = prop, data = data, selected_model = step3_model, family =family)  
  }
  
  # return selected model
  message("\n", strrep("-", getOption("width")))
  message("autoDI is limited in terms of model selection. Exercise caution when choosing your final model.")
  message(strrep("-", getOption("width")))
  return(step3_model$model)
}

AIC2 <- function(obj) {
  n <- length(na.omit(obj$y))
  p <- length(na.omit(obj$coef))
  #n <- nrow(obj$data)
  #p <- length(obj$coef)
  mu_hat <- fitted(obj)
  sigma_hat <- sqrt(sum(obj$residuals^2)/(obj$df.residual) * (n - p)/n)
  ll <- sum(dnorm(obj$y, mu_hat, sigma_hat, log = TRUE))
  np <- p + 1
  aic <- - 2*ll + 2*np
  return(aic)
}

AICc.default <- function(obj) {
  aic <- AIC2(obj)
  np <- length(obj$coef) + 1
  n <- nrow(obj$data)
  aicc <- aic + (2*np^2 + 2*np)/(n - np - 1)
  return(aicc)
}

BIC2 <- function(obj) {
  n <- length(na.omit(obj$y))
  p <- length(na.omit(obj$coef))
  #n <- nrow(obj$data)
  #p <- length(obj$coef)
  mu_hat <- fitted(obj)
  sigma_hat <- sqrt(sum(obj$residuals^2)/(obj$df.residual) * (n - p)/n)
  ll <- sum(dnorm(obj$y, mu_hat, sigma_hat, log = TRUE))
  np <- p + 1
  bic <- - 2*ll + log(n)*np
  return(bic)
}

BICc.default <- function(obj) {
  bic <- BIC2(obj)
  np <- length(obj$coef) + 1
  n <- nrow(obj$data)
  bicc <- bic + (log(n)*(np+1)*np)/(n - np - 1)
  return(bicc)
}

AICc <- function(obj) {
  UseMethod("AICc")
}

BICc <- function(obj) {
  UseMethod("BICc")
}