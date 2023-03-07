DI <- function(y, block, density, prop, treat, FG, DImodel,
               extra_formula, custom_formula,
               data, estimate_theta = FALSE, theta = 1) {
  DIcall <- match.call()
  if(theta < 0){
    stop('Please choose a positive value for theta.')
  }
  if(estimate_theta & theta != 1){
    warning('By specifying estimate_theta as TRUE, DI is overriding the specified theta value.')
  }
  
  # ensuring model tag is a string
  if(!missing(DImodel)) {
    find_input <- try(DImodel, silent = TRUE)
    if(inherits(find_input, "try-error")) {
      DImodel <- paste0(substitute(DImodel))
    } else {
      DImodel <- paste0(enquote(DImodel))[2]
    }
  }
  
  if(missing(custom_formula) & missing(y)) {
    stop("You must supply a response variable name or column index through the argument 'y'.\n")
  }
  
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
  # if model is not specified => it is a CUSTOM model
  if(missing(DImodel)) DImodel <- "CUSTOM"
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
  
  # checks for extra_formula
  if(!missing(extra_formula) & !missing(custom_formula)) {
    stop("Please provide either custom_formula or extra_formula; not the two at the same time.") 
  }
  if(missing(extra_formula)) {
    extra_formula <- 0
  }
  
  # flags if block/density and/or treat are missing
  if(missing(block)) block <- NA
  if(missing(density)) density <- NA
  treat_flag <- FALSE
  if(missing(treat)) {
    treat <- NA
    treat_flag <- TRUE
  }
  if(treat_flag & length(grep("treat", DImodel)) == 1)
    stop("you must supply the treat argument to fit model ", DImodel)
  
  # save FG names
  if(!missing(FG)) FGnames <- FG
  # NULL flag for FG variable
  if(missing(FG)) FG <- NULL
  # error message when there is no FG variable, but user specifies DImodel = "FG"
  if(DImodel == "FG" & is.null(FG)) {
    stop("The argument FG must be specified alongside DImodel = 'FG'.\n") 
  }
  if(!missing(custom_formula)) {
    if(DImodel != "CUSTOM") {
      warning("fitting custom DI model using supplied custom_formula instead of model ", DImodel,
              " (", namesub_DI(DImodel), ")")
    }
    DImodel <- "CUSTOM"
    model_fit <- DI_CUSTOM(custom_formula = custom_formula, data = data,
                           family = family, estimate_theta = estimate_theta, total = total)
  } else {
    # preparing new data object
    data_obj <- DI_data_prepare(y = y, block = block, density = density, prop = prop, treat = treat, FG = FG, data = data, theta = theta)
    newdata <- data_obj$newdata
    nSpecies <- data_obj$nSpecies
    if(nSpecies <= 2 & DImodel == "E") {
      stop("you must have > 2 species to fit model ", DImodel, " (", namesub_DI(DImodel), ")")
    }
    if(nSpecies <= 2 & DImodel == "AV") {
      stop("you must have > 2 species to fit model ", DImodel, " (", namesub_DI(DImodel), ")")
    }
    if(data_obj$P_int_flag & DImodel == "ADD") {
      stop("you must have > 3 species to fit model ", DImodel, " (", namesub_DI(DImodel), ")")
    }
    ## fitting the DI model
  
    ## adding the DI_ prefix
    model_function_name <- paste("DI_", DImodel, sep = "")
    ## adding the _treat suffix if needed
    if(!is.na(treat)) model_function_name <- paste(model_function_name, "_treat", sep = "")
    ## getting appropriate DI_ function
    model_fun <- get(model_function_name)
    ## fitting the model
    if(DImodel == "FG") {
      model_fit <- model_fun(y = data_obj$y, block = data_obj$block, density = data_obj$density, prop = data_obj$prop, FG = data_obj$FG,
                             treat = data_obj$treat, data = newdata, family = family, extra_formula = extra_formula,
                             estimate_theta = estimate_theta, nSpecies = nSpecies, total = total, FGnames = FGnames)
    } else {
      model_fit <- model_fun(y = data_obj$y, block = data_obj$block, density = data_obj$density, prop = data_obj$prop, FG = data_obj$FG,
                             treat = data_obj$treat, data = newdata, family = family, extra_formula = extra_formula,
                             estimate_theta = estimate_theta, nSpecies = nSpecies, total = total)
    }
  }
  if(!estimate_theta) {
    the_DI_model <- model_fit$model
  } else {
    the_DI_model <- model_fit$theta
  }
  message("Fitted model: ", namesub_DI(DImodel), sep = "")
  if(estimate_theta) {
    message("Theta estimate: ", round(the_DI_model$coef["theta"], 4), sep = "")  
  }
  the_DI_model$DIcall <- DIcall
  if(theta!=1 & !estimate_theta){
    the_DI_model$coefficients <- c(the_DI_model$coefficients, "theta" = theta)
    the_DI_model$df.residual <- the_DI_model$df.residual - 1
    the_DI_model$aic <- AIC2(the_DI_model)
  }
  # RV change: Adding original data to model object
  the_DI_model$original_data <- data
  #the_DI_model$aic <- AIC2(the_DI_model)
  class(the_DI_model) <- c("DI", "glm", "lm")
  return(the_DI_model)
}

DI_CUSTOM <- function(custom_formula, data, family, estimate_theta, total) {
  DIinternalcall <- match.call()
  mod_CUSTOM <- glm(formula = custom_formula, family = family, data = data)
  mod_CUSTOM$DIinternalcall <- DIinternalcall
  modlist <- list("model" = mod_CUSTOM)
  if(estimate_theta) {
    stop("theta estimation not available when custom_formula is supplied")
  }
  return(modlist)
}

DI_STR <- function(y, block, density, prop, treat, FG, data, family, estimate_theta, nSpecies, total, extra_formula = 0) {
  DIinternalcall <- match.call()
  
  if(density != "density_zero") {
    if(!inherits(extra_formula, "formula")) {
      extra_terms <- density
    } else {
      extra_terms <- paste(density, "+", paste(extra_formula)[2])
    }
  } else {
    if(!inherits(extra_formula, "formula")) {
      extra_terms <- "0"
    } else {
      extra_terms <- paste(extra_formula)[2]
    }
  }
  
  if(family %in% c("binomial","quasibinomial")) {
    y <- paste("cbind(", y, ",", total, "-", y, ")") 
  }
  if(block == "block_zero") {
    if(extra_terms == "0") {
      fmla_STR <- as.formula(paste(y, "~ 1"))
    } else {
      fmla_STR <- as.formula(paste(y, "~", extra_terms))
    }
  } else {
    if(extra_terms == "0") {
      fmla_STR <- as.formula(paste(y, "~", block))
    } else {
      fmla_STR <- as.formula(paste(y, "~", block, "+", extra_terms))
    }
  }
  mod_STR <- DI_check_and_fit(fmla = fmla_STR, y = y,
                               block = block, density = density,
                               family = family, data = data)
  mod_STR$DIinternalcall <- DIinternalcall
  return(list("model" = mod_STR,
              "theta" = mod_STR))
}

DI_STR_treat <- function(y, block, density, prop, treat, FG, data, family, estimate_theta, nSpecies, total, extra_formula = 0) {
  DIinternalcall <- match.call()
  
  if(density != "density_zero") {
    if(!inherits(extra_formula, "formula")) {
      extra_terms <- density
    } else {
      extra_terms <- paste(density, "+", paste(extra_formula)[2])
    }
  } else {
    if(!inherits(extra_formula, "formula")) {
      extra_terms <- "0"
    } else {
      extra_terms <- paste(extra_formula)[2]
    }
  }
  
  if(family %in% c("binomial","quasibinomial")) {
    y <- paste("cbind(", y, ",", total, "-", y, ")") 
  }
  if(block == "block_zero") {
    if(extra_terms == "0") {
      fmla_STR_treat <- as.formula(paste(y, "~", treat))
    } else {
      fmla_STR_treat <- as.formula(paste(y, "~", treat, "+", extra_terms))
    }
  } else {
    if(extra_terms == "0") {
      fmla_STR_treat <- as.formula(paste(y, "~", block, "+", treat))
    } else {
      fmla_STR_treat <- as.formula(paste(y, "~", block, "+", treat, "+", extra_terms))
    }
  }
  mod_STR_treat <- DI_check_and_fit(fmla = fmla_STR_treat, y = y,
                                     block = block, density = density, treat = treat,
                                     family = family, data = data)
  mod_STR_treat$DIinternalcall <- DIinternalcall
  return(list("model" = mod_STR_treat,
              "theta" = mod_STR_treat))
}

DI_ID <- function(y, block, density, prop, treat, FG, data, family, estimate_theta, nSpecies, total, extra_formula = 0) {
  DIinternalcall <- match.call()
  
  if(density != "density_zero") {
    if(!inherits(extra_formula, "formula")) {
      extra_terms <- density
    } else {
      extra_terms <- paste(density, "+", paste(extra_formula)[2])
    }
  } else {
    if(!inherits(extra_formula, "formula")) {
      extra_terms <- "0"
    } else {
      extra_terms <- paste(extra_formula)[2]
    }
  }
  
  if(family %in% c("binomial","quasibinomial")) {
    y <- paste("cbind(", y, ",", total, "-", y, ")") 
  }
  if(block == "block_zero") {
    fmla_ID <- as.formula(paste(y, "~", "0+",
                                paste(prop, collapse = "+"),
                                "+", extra_terms))
  } else {
    fmla_ID <- as.formula(paste(y, "~", "0+", paste(prop, collapse = "+"), "+", block,
                                "+", extra_terms))
  }
  # check design matrix and fit model
  mod_ID <- DI_check_and_fit(fmla = fmla_ID, y = y,
                             block = block, density = density,
                             family = family, data = data)
  mod_ID$DIinternalcall <- DIinternalcall
  return(list("model" = mod_ID,
              "theta" = mod_ID))
}

DI_ID_treat <- function(y, block, density, prop, treat, FG, data, family, estimate_theta, nSpecies, total, extra_formula = 0) {
  DIinternalcall <- match.call()
  
  if(density != "density_zero") {
    if(!inherits(extra_formula, "formula")) {
      extra_terms <- density
    } else {
      extra_terms <- paste(density, "+", paste(extra_formula)[2])
    }
  } else {
    if(!inherits(extra_formula, "formula")) {
      extra_terms <- "0"
    } else {
      extra_terms <- paste(extra_formula)[2]
    }
  }
  
  if(family %in% c("binomial","quasibinomial")) {
    y <- paste("cbind(", y, ",", total, "-", y, ")") 
  }
  if(block == "block_zero") {
    fmla_ID_treat <- as.formula(paste(y, "~", "0+", paste(prop, collapse = "+"), "+", treat,
                                      "+", extra_terms))
  } else {
    fmla_ID_treat <- as.formula(paste(y, "~", "0+", paste(prop, collapse = "+"), "+", block, "+", treat,
                                      "+", extra_terms))
  }
  # check design matrix and fit model
  mod_ID_treat <- DI_check_and_fit(fmla = fmla_ID_treat, y = y,
                                   block = block, density = density, treat = treat,
                                   family = family, data = data)
  mod_ID_treat$DIinternalcall <- DIinternalcall
  return(list("model" = mod_ID_treat,
              "theta" = mod_ID_treat))
}

DI_AV <- function(y, block, density, prop, treat, FG, data, family, estimate_theta, nSpecies, total, extra_formula = 0) {
  DIinternalcall <- match.call()
  
  if(density != "density_zero") {
    if(!inherits(extra_formula, "formula")) {
      extra_terms <- density
    } else {
      extra_terms <- paste(density, "+", paste(extra_formula)[2])
    }
  } else {
    if(!inherits(extra_formula, "formula")) {
      extra_terms <- "0"
    } else {
      extra_terms <- paste(extra_formula)[2]
    }
  }
  
  if(family %in% c("binomial","quasibinomial")) {
    y <- paste("cbind(", y, ",", total, "-", y, ")") 
  }
  if(!any(names(data) == "AV")) stop("'AV' variable not found in dataset.")
  if(block == "block_zero") {
    fmla_AV <- as.formula(paste(y, "~", "0+", 
                                paste(prop, collapse = "+"), "+AV",
                                "+", extra_terms)) 
  } else {
    fmla_AV <- as.formula(paste(y, "~", "0+",
                                paste(prop, collapse = "+"), "+AV+", block,
                                "+", extra_terms))
  }
  # check design matrix and fit model
  mod_AV <- DI_check_and_fit(fmla = fmla_AV, y = y,
                             block = block, density = density,
                             family = family, data = data)
  modlist <- list("model" = mod_AV)
  modlist$model$DIinternalcall <- DIinternalcall
  if(estimate_theta) {
    modlist$theta <- DI_theta(obj = mod_AV, DImodel = "AV",
                              nSpecies = nSpecies, family = family)
    modlist$theta$DIinternalcall <- DIinternalcall
    # RV change
    modlist$theta$DIcheck_formula <- fmla_AV
  }
  return(modlist)
}

DI_AV_treat <- function(y, block, density, prop, treat, FG, data, family, estimate_theta, nSpecies, total, extra_formula = 0) {
  DIinternalcall <- match.call()
  
  if(density != "density_zero") {
    if(!inherits(extra_formula, "formula")) {
      extra_terms <- density
    } else {
      extra_terms <- paste(density, "+", paste(extra_formula)[2])
    }
  } else {
    if(!inherits(extra_formula, "formula")) {
      extra_terms <- "0"
    } else {
      extra_terms <- paste(extra_formula)[2]
    }
  }
  
  if(family %in% c("binomial","quasibinomial")) {
    y <- paste("cbind(", y, ",", total, "-", y, ")") 
  }
  if(!any(names(data) == "AV")) stop("'AV' variable not found in dataset.")
  if(block == "block_zero") {
    fmla_AV_treat <- as.formula(paste(y, "~", "0+", paste(prop, collapse = "+"), "+AV+", treat,
                                      "+", extra_terms)) 
  } else {
    fmla_AV_treat <- as.formula(paste(y, "~", "0+", paste(prop, collapse = "+"), "+AV+", block, "+", treat,
                                      "+", extra_terms))
  }
  # check design matrix and fit model
  mod_AV_treat <- DI_check_and_fit(fmla = fmla_AV_treat, y = y,
                                   block = block, density = density, treat = treat,
                                   family = family, data = data)
  modlist <- list("model" = mod_AV_treat)
  modlist$model$DIinternalcall <- DIinternalcall
  if(estimate_theta) {
    modlist$theta <- DI_theta(obj = mod_AV_treat, DImodel = "AV",
                              nSpecies = nSpecies, family = family)
    modlist$theta$DIinternalcall <- DIinternalcall
    # RV change
    modlist$theta$DIcheck_formula <- fmla_AV_treat
  }
  return(modlist)
}

DI_E <- function(y, block, density, prop, treat, data, FG, family, estimate_theta, nSpecies, total, extra_formula = 0) {
  DIinternalcall <- match.call()
  
  if(density != "density_zero") {
    if(!inherits(extra_formula, "formula")) {
      extra_terms <- density
    } else {
      extra_terms <- paste(density, "+", paste(extra_formula)[2])
    }
  } else {
    if(!inherits(extra_formula, "formula")) {
      extra_terms <- "0"
    } else {
      extra_terms <- paste(extra_formula)[2]
    }
  }
  
  if(family %in% c("binomial","quasibinomial")) {
    y <- paste("cbind(", y, ",", total, "-", y, ")") 
  }
  if(!any(names(data) == "E")) stop("'E' variable not found in dataset.")
  if(block == "block_zero") {
    fmla_E <- as.formula(paste(y, "~", "0+", paste(prop, collapse = "+"), "+E",
                               "+", extra_terms)) 
  } else {
    fmla_E <- as.formula(paste(y, "~", "0+", paste(prop, collapse = "+"), "+E+", block,
                               "+", extra_terms))
  }
  # check design matrix and fit model
  mod_E <- DI_check_and_fit(fmla = fmla_E, y = y,
                            block = block, density = density,
                            family = family, data = data)
  modlist <- list("model" = mod_E)
  modlist$model$DIinternalcall <- DIinternalcall
  if(estimate_theta) {
    modlist$theta <- DI_theta(obj = mod_E, DImodel = "E",
                              nSpecies = nSpecies, family = family)
    modlist$theta$DIinternalcall <- DIinternalcall
    # RV change
    modlist$theta$DIcheck_formula <- fmla_E
  }
  return(modlist)
}

DI_E_treat <- function(y, block, density, prop, treat, data, FG, family, estimate_theta, nSpecies, total, extra_formula = 0) {
  DIinternalcall <- match.call()
  
  if(density != "density_zero") {
    if(!inherits(extra_formula, "formula")) {
      extra_terms <- density
    } else {
      extra_terms <- paste(density, "+", paste(extra_formula)[2])
    }
  } else {
    if(!inherits(extra_formula, "formula")) {
      extra_terms <- "0"
    } else {
      extra_terms <- paste(extra_formula)[2]
    }
  }
  
  if(family %in% c("binomial","quasibinomial")) {
    y <- paste("cbind(", y, ",", total, "-", y, ")") 
  }
  if(!any(names(data) == "E")) stop("'E' variable not found in dataset.")
  if(block == "block_zero") {
    fmla_E_treat <- as.formula(paste(y, "~", "0+", paste(prop, collapse = "+"), "+E+", treat,
                                     "+", extra_terms)) 
  } else {
    fmla_E_treat <- as.formula(paste(y, "~", "0+", paste(prop, collapse = "+"), "+E+", block, "+", treat,
                                     "+", extra_terms))
  }
  # check design matrix and fit model
  mod_E_treat <- DI_check_and_fit(fmla = fmla_E_treat, y = y,
                                  block = block, density = density, treat = treat,
                                  family = family, data = data)
  modlist <- list("model" = mod_E_treat)
  modlist$model$DIinternalcall <- DIinternalcall
  if(estimate_theta) {
    modlist$theta <- DI_theta(obj = mod_E_treat, DImodel = "E",
                              nSpecies = nSpecies, family = family)
    modlist$theta$DIinternalcall <- DIinternalcall 
    # RV change
    modlist$theta$DIcheck_formula <- fmla_E_treat
  }
  return(modlist)
}

DI_FG <- function(y, block, density, prop, treat, FG, data, family, estimate_theta, nSpecies, total, FGnames, extra_formula = 0) {
  DIinternalcall <- match.call()
  
  if(density != "density_zero") {
    if(!inherits(extra_formula, "formula")) {
      extra_terms <- density
    } else {
      extra_terms <- paste(density, "+", paste(extra_formula)[2])
    }
  } else {
    if(!inherits(extra_formula, "formula")) {
      extra_terms <- "0"
    } else {
      extra_terms <- paste(extra_formula)[2]
    }
  }
  
  if(family %in% c("binomial","quasibinomial")) {
    y <- paste("cbind(", y, ",", total, "-", y, ")") 
  }
  if(is.null(FG)) stop("FG matrix not found")
  #FG_ <- FG
  if(block == "block_zero") {
    fmla_FG <- as.formula(paste(y, "~", "0+", paste(prop, collapse = "+"), "+FG_",
                                "+", extra_terms))
  } else {
    fmla_FG <- as.formula(paste(y, "~", "0+", paste(prop, collapse = "+"), "+FG_+", block,
                                "+", extra_terms))
  }
  # check design matrix and fit model
  mod_FG <- DI_check_and_fit(fmla = fmla_FG, y = y, FG = FG,
                             block = block, density = density,
                             family = family, data = data)
  modlist <- list("model" = mod_FG)
  modlist$model$DIinternalcall <- DIinternalcall
  if(estimate_theta) {
    modlist$theta <- DI_theta(obj = mod_FG, DImodel = "FG", FGnames = FGnames,
                              nSpecies = nSpecies, family = family)
    modlist$theta$DIinternalcall <- DIinternalcall
    # RV change
    modlist$theta$DIcheck_formula <- fmla_FG
  }
  return(modlist)
}

DI_FG_treat <- function(y, block, density, prop, treat, FG, data, family, estimate_theta, nSpecies, total, FGnames, extra_formula = 0) {
  DIinternalcall <- match.call()
  
  if(density != "density_zero") {
    if(!inherits(extra_formula, "formula")) {
      extra_terms <- density
    } else {
      extra_terms <- paste(density, "+", paste(extra_formula)[2])
    }
  } else {
    if(!inherits(extra_formula, "formula")) {
      extra_terms <- "0"
    } else {
      extra_terms <- paste(extra_formula)[2]
    }
  }
  
  if(family %in% c("binomial","quasibinomial")) {
    y <- paste("cbind(", y, ",", total, "-", y, ")") 
  }
  if(is.null(FG)) stop("FG matrix not found")
  #FG_ <- FG
  if(block == "block_zero") {
    fmla_FG_treat <- as.formula(paste(y, "~", "0+", paste(prop, collapse = "+"), "+FG_+", treat,
                                      "+", extra_terms))
  } else {
    fmla_FG_treat <- as.formula(paste(y, "~", "0+", paste(prop, collapse = "+"), "+FG_+", block, "+", treat,
                                      "+", extra_terms))
  }
  # check design matrix and fit model
  mod_FG_treat <- DI_check_and_fit(fmla = fmla_FG_treat, y = y, FG = FG,
                                   block = block, density = density, treat = treat,
                                   family = family, data = data)
  modlist <- list("model" = mod_FG_treat)
  modlist$model$DIinternalcall <- DIinternalcall
  if(estimate_theta) {
    modlist$theta <- DI_theta(obj = mod_FG_treat, DImodel = "FG", FGnames = FGnames,
                              nSpecies = nSpecies, family = family)
    modlist$theta$DIinternalcall <- DIinternalcall
    # RV change
    modlist$theta$DIcheck_formula <- fmla_FG_treat
  }
  return(modlist)
}

DI_ADD <- function(y, block, density, prop, treat, FG, data, family, estimate_theta, nSpecies, total, extra_formula = 0) {
  DIinternalcall <- match.call()
  
  if(density != "density_zero") {
    if(!inherits(extra_formula, "formula")) {
      extra_terms <- density
    } else {
      extra_terms <- paste(density, "+", paste(extra_formula)[2])
    }
  } else {
    if(!inherits(extra_formula, "formula")) {
      extra_terms <- "0"
    } else {
      extra_terms <- paste(extra_formula)[2]
    }
  }
  
  if(family %in% c("binomial","quasibinomial")) {
    y <- paste("cbind(", y, ",", total, "-", y, ")") 
  }
  nam <- names(data)
  vec_grep <- Vectorize(grep, "pattern")
  if(length(vec_grep(paste0(prop, "_add"), nam)) != length(prop)) 
    stop("P_int variables not found in dataset.")
  if(block == "block_zero") {
    fmla_ADD <- as.formula(paste(y, "~", "0+", 
                                    paste(prop, collapse = "+"), "+",
                                    paste0(prop, "_add", collapse = "+"),
                                 "+", extra_terms))
  } else {
    fmla_ADD <- as.formula(paste(y, "~", "0+",
                                    paste(prop, collapse = "+"), "+",
                                    paste0(prop, "_add", collapse = "+"), "+", block,
                                 "+", extra_terms))
  }
  # check design matrix and fit model
  mod_ADD <- DI_check_and_fit(fmla = fmla_ADD, y = y,
                              block = block, density = density,
                              family = family, data = data)
  modlist <- list("model" = mod_ADD)
  modlist$model$DIinternalcall <- DIinternalcall
  if(estimate_theta) {
    modlist$theta <- DI_theta(obj = mod_ADD, DImodel = "ADD",
                              nSpecies = nSpecies, family = family)
    modlist$theta$DIinternalcall <- DIinternalcall
    # RV change
    modlist$theta$DIcheck_formula <- fmla_ADD
  }
  return(modlist)
}

DI_ADD_treat <- function(y, block, density, prop, treat, FG, data, family, estimate_theta, nSpecies, total, extra_formula = 0) {
  DIinternalcall <- match.call()
  
  if(density != "density_zero") {
    if(!inherits(extra_formula, "formula")) {
      extra_terms <- density
    } else {
      extra_terms <- paste(density, "+", paste(extra_formula)[2])
    }
  } else {
    if(!inherits(extra_formula, "formula")) {
      extra_terms <- "0"
    } else {
      extra_terms <- paste(extra_formula)[2]
    }
  }
  
  if(family %in% c("binomial","quasibinomial")) {
    y <- paste("cbind(", y, ",", total, "-", y, ")") 
  }
  nam <- names(data)
  vec_grep <- Vectorize(grep, "pattern")
  if(length(vec_grep(paste0(prop, "_add"), nam)) != length(prop)) 
    stop("P_int variables not found in dataset.")
  if(block == "block_zero") {
    fmla_ADD_treat <- as.formula(paste(y, "~", "0+",
                                     paste(prop, collapse = "+"), "+",
                                     paste0(prop, "_add", collapse = "+"), "+", treat,
                                     "+", extra_terms))
  } else {
    fmla_ADD_treat <- as.formula(paste(y, "~", "0+", 
                                     paste(prop, collapse = "+"), "+",
                                     paste0(prop, "_add", collapse = "+"), "+", block, "+", treat,
                                     "+", extra_terms))
  }
  # check design matrix and fit model
  mod_ADD_treat <- DI_check_and_fit(fmla = fmla_ADD_treat, y = y,
                                    block = block, density = density, treat = treat,
                                    family = family, data = data)
  modlist <- list("model" = mod_ADD_treat)
  modlist$model$DIinternalcall <- DIinternalcall
  if(estimate_theta) {
    modlist$theta <- DI_theta(obj = mod_ADD_treat, DImodel = "ADD",
                              nSpecies = nSpecies, family = family)
    modlist$theta$DIinternalcall <- DIinternalcall
    # RV change
    modlist$theta$DIcheck_formula <- fmla_ADD_treat
  }
  return(modlist)
}

DI_FULL <- function(y, block, density, prop, treat, FG, data, family, estimate_theta, nSpecies, total, extra_formula = 0) {
  DIinternalcall <- match.call()
  
  # Getting the FULL variables
  #print(DI_data_fullpairwise(prop = prop, data = data[1, prop], theta = 1))
  FULL <- DI_data_fullpairwise(prop = prop, data = data[1, prop], theta = 1)
  full_names <- paste0('`', names(FULL), '`')
  
  if(density != "density_zero") {
    if(!inherits(extra_formula, "formula")) {
      extra_terms <- density
    } else {
      extra_terms <- paste(density, "+", paste(extra_formula)[2])
    }
  } else {
    if(!inherits(extra_formula, "formula")) {
      extra_terms <- "0"
    } else {
      extra_terms <- paste(extra_formula)[2]
    }
  }
  
  if(family %in% c("binomial","quasibinomial")) {
    y <- paste("cbind(", y, ",", total, "-", y, ")") 
  }
  if(block == "block_zero") {
    fmla_FULL <- as.formula(paste(y, "~", "0+", 
                                  paste(prop, collapse = "+"), "+",
                                  paste0(full_names, collapse = "+"),
                                  "+", extra_terms))
  } else {
    fmla_FULL <- as.formula(paste(y, "~", "0+",
                                  paste(prop, collapse = "+"), "+",
                                  paste0(full_names, collapse = "+"), "+", block,
                                  "+", extra_terms))
  }
  # check design matrix and fit model
  mod_FULL <- DI_check_and_fit(fmla = fmla_FULL, y = y,
                               block = block, density = density,
                               family = family, data = data)
  modlist <- list("model" = mod_FULL)
  modlist$model$DIinternalcall <- DIinternalcall
  if(estimate_theta) {
    modlist$theta <- DI_theta(obj = mod_FULL, DImodel = "FULL", prop = prop,
                              nSpecies = nSpecies, family = family)
    modlist$theta$DIinternalcall <- DIinternalcall
    # RV change
    modlist$theta$DIcheck_formula <- fmla_FULL
  }
  return(modlist)
}

DI_FULL_treat <- function(y, block, density, prop, treat, FG, data, family, estimate_theta, nSpecies, total, extra_formula = 0) {
  DIinternalcall <- match.call()
  
  # Getting the FULL variables
  FULL <- DI_data_fullpairwise(prop = prop, data = data[1, prop], theta = 1)
  full_names <- paste0('`', names(FULL), '`')
  
  if(density != "density_zero") {
    if(!inherits(extra_formula, "formula")) {
      extra_terms <- density
    } else {
      extra_terms <- paste(density, "+", paste(extra_formula)[2])
    }
  } else {
    if(!inherits(extra_formula, "formula")) {
      extra_terms <- "0"
    } else {
      extra_terms <- paste(extra_formula)[2]
    }
  }
  
  if(family %in% c("binomial","quasibinomial")) {
    y <- paste("cbind(", y, ",", total, "-", y, ")") 
  }
  if(block == "block_zero") {
    fmla_FULL_treat <- as.formula(paste(y, "~", "0+",
                                        paste(prop, collapse = "+"), "+",
                                        paste0(full_names, collapse = "+"), "+", treat,
                                        "+", extra_terms))
  } else {
    fmla_FULL_treat <- as.formula(paste(y, "~", "0+", 
                                        paste(prop, collapse = "+"), "+",
                                        paste0(full_names, collapse = "+"), "+", block, "+", treat,
                                        "+", extra_terms))
  }
  # check design matrix and fit model
  mod_FULL_treat <- DI_check_and_fit(fmla = fmla_FULL_treat, y = y,
                                     block = block, density = density, treat = treat,
                                     family = family, data = data)
  modlist <- list("model" = mod_FULL_treat)
  modlist$model$DIinternalcall <- DIinternalcall
  if(estimate_theta) {
    modlist$theta <- DI_theta(obj = mod_FULL_treat, DImodel = "FULL", prop = prop,
                              nSpecies = nSpecies, family = family)
    modlist$theta$DIinternalcall <- DIinternalcall
    # RV change
    modlist$theta$DIcheck_formula <- fmla_FULL_treat
  }
  return(modlist)
}

#################################
## S3 methods for DI class ##
#################################

AIC.DI <- function(object, ...) {
  AIC2(object)
}

BIC.DI <- function(object, ...) {
  BIC2(object)
}

AICc.DI <- function(obj) {
  AICc.default(obj)
}

BICc.DI <- function(obj) {
  BICc.default(obj)
}

logLik.DI <- function(object, ...) {
  if (!missing(...)) 
    warning("extra arguments discarded")
  fam <- family(object)$family
  p <- object$df.null - object$df.residual
  if (fam %in% c("gaussian", "Gamma", "inverse.gaussian")) 
    p <- p + 1
  val <- p - object$aic/2
  attr(val, "nobs") <- sum(!is.na(object$residuals))
  attr(val, "df") <- p
  class(val) <- "logLik"
  val
}

extract <- function(obj, what = c("FULL","ADD","FG","AV","E")) {
  UseMethod("extract")
}

extract.DI <- function(obj, what = c("FULL","ADD","FG","AV","E")) {
  all_vars <- names(obj$data)
  ret <- list()
  if("E" %in% what) {
    all_vars_end <- substr(all_vars, start = nchar(all_vars) - 1, stop = nchar(all_vars))
    E_vars <- grep("E", all_vars_end)
    ret$E <- obj$data[,E_vars]
  }
  if("AV" %in% what) {
    all_vars_end <- substr(all_vars, start = nchar(all_vars) - 2, stop = nchar(all_vars))
    AV_vars <- grep("AV", all_vars_end)
    ret$AV <- obj$data[,AV_vars]
  }
  if("FG" %in% what) {
    FG_vars <- grep("FG", all_vars)
    ret$FG <- obj$data[,FG_vars]
  }
  if("ADD" %in% what) {
    all_vars_end <- substr(all_vars, start = nchar(all_vars) - 3, stop = nchar(all_vars))
    add_vars <- grep("_add", all_vars_end)
    ret$ADD <- obj$data[,add_vars]
  }
  if("FULL" %in% what) {
    all_vars <- names(obj$data)
    pairwise_vars <- grep(":", all_vars)
    ret$pairwise <- obj$data[,pairwise_vars]
  }
  return(ret)
}

extract.autoDI <- function(obj, what = c("FULL","ADD","FG","AV","E")) {
  new_obj <- obj$selected_model_obj
  extract.DI(obj = new_obj, what = what)
}