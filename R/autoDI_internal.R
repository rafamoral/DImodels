namesub_autoDI <- Vectorize(function(name) {
  thename <- switch(name,
                    "STR_model" = "Structural 'STR' DImodel",
                    "ID_model" = "Species identity 'ID' DImodel",
                    "FULL_model" = "Separate pairwise interactions 'FULL' DImodel",
                    "AV_model" = "Average interactions 'AV' DImodel",
                    "E_model" = "Evenness 'E' DImodel",
                    "ADD_model" = 
                      "Additive species contributions to interactions 'ADD' DImodel",
                    "FG_model" = "Functional group effects 'FG' DImodel",
                    "STR_model_treat" = "Structural 'STR' DImodel with treatment",
                    "ID_model_treat" = "Species identity 'ID' DImodel with treatment",
                    "FULL_model_treat" = "Separate pairwise interactions 'FULL' DImodel with treatment",
                    "AV_model_treat" = "Average interactions 'AV' DImodel with treatment",
                    "E_model_treat" = "Evenness 'E' DImodel with treatment",
                    "ADD_model_treat" = 
                      "Additive species contributions to interactions 'ADD' DImodel with treatment",
                    "FG_model_treat" = "Functional group effects 'FG' DImodel with treatment",
                    "FULL_model_theta" = "Separate pairwise interactions 'FULL' DImodel, estimating theta",
                    "AV_model_theta" = "Average interactions 'AV' DImodel, estimating theta",
                    "E_model_theta" = "Evenness 'E' DImodel, estimating theta",
                    "ADD_model_theta" = 
                      "Additive species contributions to interactions 'ADD' DImodel, estimating theta",
                    "FG_model_theta" = "Functional group effects 'FG' DImodel, estimating theta",
                    "FULL_model_treat_theta" = "Separate pairwise interactions 'FULL' DImodel with treatment, estimating theta",
                    "AV_model_treat_theta" = "Average interactions 'AV' DImodel with treatment, estimating theta",
                    "E_model_treat_theta" = "Evenness 'E' DImodel with treatment, estimating theta",
                    "ADD_model_treat_theta" = 
                      "Additive species contributions to interactions 'ADD' DImodel with treatment, estimating theta",
                    "FG_model_treat_theta" = "Functional group effects 'FG' DImodel with treatment, estimating theta",
                    stop("not yet implemented"))
  return(thename)
}, "name")

reftest_autoDI <- function(model_to_compare, ref_model, family) {
  if(family %in% c("poisson","binomial")) {
    ref_test <- "Chisq"
  } else {
    ref_test <- "F"
  }
  anovas <- anova(model_to_compare, ref_model, test = ref_test)
  ## formatting the output
  anovas_format <- as.data.frame(anovas)
  anovas_format <- round(anovas_format, 4)
  anovas_format[is.na(anovas_format)] <- ""
  names(anovas_format)[c(2,4)] <- c("Resid. SSq","SSq")
  anovas_format$"Resid. MSq" <- round(anovas_format[,2]/anovas_format[,1], 4)
  anovas_format <- anovas_format[,c(1,2,7,3:6)]
  model_tokens <- c("Selected","Reference")
  anovas_format$model <- model_tokens
  anovas_format <- anovas_format[,c(8,1:7)]
  row.names(anovas_format) <- paste("DI Model", row.names(anovas_format))
  anovas_format[,8][anovas_format[,8] == 0] <- "<0.0001"
  message("\n")
  old <- options()
  options(scipen = 999)
  on.exit(options(old))
  print(anovas_format)
  #options(scipen = 0)
  #message("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
}

test_autoDI <- function(model_list, family, treat) {
  if(family %in% c("poisson","binomial")) {
    message("Selection using X2 tests", "\n")
    Test <- "Chisq"
  } else {
    message("Selection using F tests", "\n")
    Test <- "F"
  }
  
  model_names <- names(model_list)
  treat_flag <- grep("treatment", namesub_autoDI(model_names))
  theta_flag <- grep("theta", namesub_autoDI(model_names))
  
  treat_output <- rep("none", length(model_names))
  treat_output[treat_flag] <- paste("'", treat, "'", sep = "")
  theta_output <- rep(FALSE, length(model_names))
  theta_output[theta_flag] <- TRUE
  
  model_tokens <- gsub("_", "", gsub("[:a-z:]", "", names(model_list)))
  names(model_list) <- NULL
  #anovas <- do.call(anova, model_list) ## if using lm objects
  anovas <- eval(parse(text = paste("anova(",
                                    paste("model_list[[", 1:length(model_list), "]]",
                                          sep = "", collapse = ","),
                                    ",test ='", Test, "')", sep = "")
                                   ))
  if(family %in% c("poisson","binomial")) {
    p_values <- anovas$"Pr(>Chi)"
  } else {
    p_values <- anovas$"Pr(>F)"
  }
  p_less <- which(p_values < .05)
  p_value_selected <- ifelse(length(p_less) == 0, 1, max(p_less))
  selected <- model_names[p_value_selected]
  ## formatting the output
  anovas_format <- as.data.frame(anovas)
  anovas_format <- round(anovas_format, 4)
  anovas_format[is.na(anovas_format)] <- ""
  names(anovas_format)[c(2,4)] <- c("Resid. SSq","SSq")
  anovas_format$"Resid. MSq" <- round(anovas_format[,2]/anovas_format[,1], 4)
  anovas_format <- anovas_format[,c(1,2,7,3:6)]
  anovas_format$model <- model_tokens
  anovas_format$treat <- treat_output
  anovas_format$theta <- theta_output
  anovas_format <- anovas_format[,c(8:10,1:7)]
  names(anovas_format)[1] <- "DI_model"
  names(anovas_format)[3] <- "estimate_theta"
  row.names(anovas_format) <- paste("DI Model", row.names(anovas_format))
  anovas_format[,10][anovas_format[,10] == 0] <- "<0.0001"
  desc_table <- data.frame("Description" = namesub_autoDI(model_names))
  row.names(desc_table) <- paste("DI Model", 1:nrow(anovas_format))
  print(desc_table, right = FALSE)
  message("\n")
  old <- options()
  options(scipen = 999)
  on.exit(options(old))
  print(anovas_format)
  #options(scipen = 0)
  #message("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  return(selected)
}

  AICsel_autoDI <- function(model_list, mAIC, treat) {
  message("Selection by AIC\nWarning: DI Model with the lowest AIC will be selected, even if the difference is very small.\nPlease inspect other models to see differences in AIC.", "\n\n", sep = "")
  
  model_names <- names(model_list)
  treat_flag <- grep("treatment", namesub_autoDI(model_names))
  theta_flag <- grep("theta", namesub_autoDI(model_names))
  
  treat_output <- rep("none", length(model_names))
  treat_output[treat_flag] <- paste("'", treat, "'", sep = "")
  theta_output <- rep(FALSE, length(model_names))
  theta_output[theta_flag] <- TRUE
  
  model_tokens <- gsub("_", "", gsub("[:a-z:]", "", model_names))
  model_descriptions <- namesub_autoDI(model_names)
  the_table <- data.frame("AIC" = mAIC,
                          "DI_Model" = model_tokens,
                          "treat" = treat_output,
                          "theta" = theta_output,
                          "Description" = model_descriptions,
                          row.names = NULL)
  print(the_table, right = FALSE)
  selected <- names(model_list)[which.min(mAIC)]
  return(selected)
}

AICcsel_autoDI <- function(model_list, mAICc, treat) {
  message("Selection by AICc\nWarning: DI Model with the lowest AICc will be selected, even if the difference is very small.\nPlease inspect other models to see differences in AICc.", "\n\n", sep = "")
  
  model_names <- names(model_list)
  treat_flag <- grep("treatment", namesub_autoDI(model_names))
  theta_flag <- grep("theta", namesub_autoDI(model_names))
  
  treat_output <- rep("none", length(model_names))
  treat_output[treat_flag] <- paste("'", treat, "'", sep = "")
  theta_output <- rep(FALSE, length(model_names))
  theta_output[theta_flag] <- TRUE
  
  model_tokens <- gsub("_", "", gsub("[:a-z:]", "", model_names))
  model_descriptions <- namesub_autoDI(model_names)
  the_table <- data.frame("AICc" = mAICc,
                          "DI_Model" = model_tokens,
                          "treat" = treat_output,
                          "theta" = theta_output,
                          "Description" = model_descriptions,
                          row.names = NULL)
  print(the_table, right = FALSE)
  selected <- names(model_list)[which.min(mAICc)]
  return(selected)
}

BICsel_autoDI <- function(model_list, mBIC, treat) {
  message("Selection by BIC\nWarning: DI Model with the lowest BIC will be selected, even if the difference is very small.\nPlease inspect other models to see differences in BIC.", "\n\n", sep = "")
  
  model_names <- names(model_list)
  treat_flag <- grep("treatment", namesub_autoDI(model_names))
  theta_flag <- grep("theta", namesub_autoDI(model_names))
  
  treat_output <- rep("none", length(model_names))
  treat_output[treat_flag] <- paste("'", treat, "'", sep = "")
  theta_output <- rep(FALSE, length(model_names))
  theta_output[theta_flag] <- TRUE
  
  model_tokens <- gsub("_", "", gsub("[:a-z:]", "", model_names))
  model_descriptions <- namesub_autoDI(model_names)
  the_table <- data.frame("BIC" = mBIC,
                          "DI_Model" = model_tokens,
                          "treat" = treat_output,
                          "theta" = theta_output,
                          "Description" = model_descriptions,
                          row.names = NULL)
  print(the_table, right = FALSE)
  selected <- names(model_list)[which.min(mBIC)]
  return(selected)
}

BICcsel_autoDI <- function(model_list, mBICc, treat) {
  message("Selection by BICc\nWarning: DI Model with the lowest BICc will be selected, even if the difference is very small.\nPlease inspect other models to see differences in BICc.", "\n\n", sep = "")
  
  model_names <- names(model_list)
  treat_flag <- grep("treatment", namesub_autoDI(model_names))
  theta_flag <- grep("theta", namesub_autoDI(model_names))
  
  treat_output <- rep("none", length(model_names))
  treat_output[treat_flag] <- paste("'", treat, "'", sep = "")
  theta_output <- rep(FALSE, length(model_names))
  theta_output[theta_flag] <- TRUE
  
  model_tokens <- gsub("_", "", gsub("[:a-z:]", "", model_names))
  model_descriptions <- namesub_autoDI(model_names)
  the_table <- data.frame("BICc" = mBICc,
                          "DI_Model" = model_tokens,
                          "treat" = treat_output,
                          "theta" = theta_output,
                          "Description" = model_descriptions,
                          row.names = NULL)
  print(the_table, right = FALSE)
  selected <- names(model_list)[which.min(mBICc)]
  return(selected)
}

DI_matrix_check <- function(model_matrix) {
  n_parms <- ncol(model_matrix)
  matrix_rank <- qr(model_matrix)$rank
  if(matrix_rank == n_parms) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

autoDI_step0 <- function(y, block, density, treat, family, data) {
  if(is.na(block) & is.na(density) & is.na(treat)) {
    return(invisible())
  }
  
  message("\n", strrep("-", getOption("width")))
  message("\nSequential analysis: Investigating only non-diversity experimental design structures\n")
  fmla1 <- paste(y, "~", 1)
  fit1 <- glm(fmla1, family = family, data = data)
  
  if(!is.na(block) & !is.na(density) & !is.na(treat)) {
    fmla2 <- paste(y, "~", block)
    fmla3 <- paste(y, "~", block, "+", density)
    fmla4 <- paste(y, "~", block, "+", density, "+", treat)
    fit2 <- glm(fmla2, family = family, data = data)
    fit3 <- glm(fmla3, family = family, data = data)
    fit4 <- glm(fmla4, family = family, data = data)
    model_list <- list(fit1, fit2, fit3, fit4)
    model_tokens <- c("Intercept only","block","block + density","block + density + treat")
  }
  
  if(!is.na(block) & !is.na(density) & is.na(treat)) {
    fmla2 <- paste(y, "~", block)
    fmla3 <- paste(y, "~", block, "+", density)
    fit2 <- glm(fmla2, family = family, data = data)
    fit3 <- glm(fmla3, family = family, data = data)
    model_list <- list(fit1, fit2, fit3)
    model_tokens <- c("Intercept only","block","block + density")
  }
  
  if(!is.na(block) & is.na(density) & !is.na(treat)) {
    fmla2 <- paste(y, "~", block)
    fmla3 <- paste(y, "~", block, "+", treat)
    fit2 <- glm(fmla2, family = family, data = data)
    fit3 <- glm(fmla3, family = family, data = data)
    model_list <- list(fit1, fit2, fit3)
    model_tokens <- c("Intercept only","block","block + treat")
  }
  
  if(is.na(block) & !is.na(density) & !is.na(treat)) {
    fmla2 <- paste(y, "~", density)
    fmla3 <- paste(y, "~", density, "+", treat)
    fit2 <- glm(fmla2, family = family, data = data)
    fit3 <- glm(fmla3, family = family, data = data)
    model_list <- list(fit1, fit2, fit3)
    model_tokens <- c("Intercept only","density","density + treat")
  }
  
  if(!is.na(block) & is.na(density) & is.na(treat)) {
    fmla2 <- paste(y, "~", block)
    fit2 <- glm(fmla2, family = family, data = data)
    model_list <- list(fit1, fit2)
    model_tokens <- c("Intercept only","block")
  }
  
  if(is.na(block) & !is.na(density) & is.na(treat)) {
    fmla2 <- paste(y, "~", density)
    fit2 <- glm(fmla2, family = family, data = data)
    model_list <- list(fit1, fit2)
    model_tokens <- c("Intercept only","density")
  }
  
  if(is.na(block) & is.na(density) & !is.na(treat)) {
    fmla2 <- paste(y, "~", treat)
    fit2 <- glm(fmla2, family = family, data = data)
    model_list <- list(fit1, fit2)
    model_tokens <- c("Intercept only","treat")
  }
  
  if(family %in% c("poisson","binomial")) {
    Test <- "Chisq"
  } else {
    Test <- "F"
  }
  
  anovas <- eval(parse(text = paste("anova(",
                                    paste("model_list[[", 1:length(model_list), "]]",
                                          sep = "", collapse = ","),
                                    ",test ='", Test, "')", sep = "")
  ))
  
  ## formatting the output
  anovas_format <- as.data.frame(anovas)
  anovas_format <- round(anovas_format, 4)
  anovas_format[is.na(anovas_format)] <- ""
  names(anovas_format)[c(2,4)] <- c("Resid. SSq","SSq")
  anovas_format$"Resid. MSq" <- round(anovas_format[,2]/anovas_format[,1], 4)
  anovas_format <- anovas_format[,c(1,2,7,3:6)]
  anovas_format$model <- model_tokens
  anovas_format <- anovas_format[,c(8,1:7)]
  anovas_format[,8][anovas_format[,8] == 0] <- "<0.0001"
  message("\n")
  old <- options()
  options(scipen = 999)
  on.exit(options(old))
  print(anovas_format)
  #options(scipen = 0)
}

autoDI_step1 <- function(y, block, density, prop, treat, family, data, selection) {
  
  model_tag <- ifelse(length(prop) == 2, "FULL", "AV")
  
  if(is.na(treat)) {
      model_name <- paste0(model_tag, "_model")
      fit_theta <- suppressWarnings(suppressMessages(DI(y = y, block = block, density = density, prop = prop,
                                       DImodel = model_tag, data = data, estimate_theta = TRUE)))
      fit_notheta <- suppressWarnings(suppressMessages(DI(y = y, block = block, density = density, prop = prop,
                                         DImodel = model_tag, data = data, estimate_theta = FALSE)))
      fit_theta$DIcall <- call("DI", y = y, block = block, density = density, prop = prop,
                               DImodel = model_tag, data = data, estimate_theta = TRUE)
      fit_notheta$DIcall <- call("DI", y = y, block = block, density = density, prop = prop,
                                 DImodel = model_tag, data = data, estimate_theta = FALSE)
    } else {
      model_name <- paste0(model_tag, "_model_treat")
      fit_theta <- suppressWarnings(suppressMessages(DI(y = y, block = block, density = density, prop = prop,
                                       treat = treat, DImodel = model_tag, data = data, estimate_theta = TRUE)))
      fit_notheta <- suppressWarnings(suppressMessages(DI(y = y, block = block, density = density, prop = prop,
                                         treat = treat, DImodel = model_tag, data = data, estimate_theta = FALSE)))
      fit_theta$DIcall <- call("DI", y = y, block = block, density = density, prop = prop,
                               treat = treat, DImodel = model_tag, data = data, estimate_theta = TRUE)
      fit_notheta$DIcall <- call("DI", y = y, block = block, density = density, prop = prop,
                                 treat = treat, DImodel = model_tag, data = data, estimate_theta = FALSE)
    }
  
  message("\n", strrep("-", getOption("width")))
  message("Step 1: Investigating whether theta is equal to 1 or not for the ", model_tag, " model, including all available structures")
  th <- fit_theta$coef["theta"]
  message("\nTheta estimate: ", round(th, 4), "\n", sep = "")
  final_model_list <- list(fit_notheta,
                           fit_theta)
  names(final_model_list) <- c(model_name,
                               paste(model_name, "_theta", sep = ""))
  
  mAIC_final <- sapply(final_model_list, AIC2)
  mAICc_final <- sapply(final_model_list, AICc)
  mBIC_final <- sapply(final_model_list, BIC2)
  mBICc_final <- sapply(final_model_list, BICc)
  selected_model_final <- switch(selection,
                                 Ftest = test_autoDI(model_list = final_model_list, family = family, treat = treat),
                                 AIC = AICsel_autoDI(model_list = final_model_list, mAIC = mAIC_final, treat = treat),
                                 AICc = AICcsel_autoDI(model_list = final_model_list, mAICc = mAICc_final, treat = treat),
                                 BIC = BICsel_autoDI(model_list = final_model_list, mBIC = mBIC_final, treat = treat),
                                 BICc = BICcsel_autoDI(model_list = final_model_list, mBICc = mBICc_final, treat = treat)
  )
  
  theta_flag <- length(grep("theta", selected_model_final)) == 1
  
  if(theta_flag) conclusion <- "" else conclusion <- "not "
  
  message("\nThe test concludes that theta is ", conclusion, "significantly different from 1.")  
  
  if(theta_flag) fit_final <- fit_theta else fit_final <- fit_notheta
  
  return(list("model" = fit_final,
              "model_name" = selected_model_final,
              "theta_flag" = theta_flag))
}

autoDI_step2 <- function(y, block, density, 
                         prop, treat, FG, 
                         family, data, theta_flag,
                         selection, selected_model){
  
  message("\n", strrep("-", getOption("width")))
  message("Step 2: Investigating the interactions\n")
  if(!is.na(block) & is.na(density)) message("All models include block\n")
  if(is.na(block) & !is.na(density)) message("All models include density\n")
  if(!is.na(block) & !is.na(density)) message("All models include block and density\n")
  
  AV_flag <- TRUE
  ADD_flag <- TRUE
  FG_flag <- TRUE
  FULL_flag <- TRUE
  
  if(length(prop) == 2){
    AV_flag <- FALSE
    ADD_flag <- FALSE
    FG_flag <- FALSE
  }
  if(length(prop) ==3){
    ADD_flag <- FALSE
  }
  if(unique(is.null(FG))){
    FG_flag <-  FALSE
  }
  
  ID_model <- suppressWarnings(suppressMessages(DI(y = y, block = block, density = density, 
                 prop = prop, treat = treat, 
                 data = data, DImodel = 'ID')))
  STR_model <- suppressWarnings(suppressMessages(DI(y = y, block = block, density = density, 
                  prop = prop, treat = treat, 
                  data = data, DImodel = 'STR')))
  model_list <- list('STR_model' = STR_model, 'ID_model' = ID_model)
  
  if(theta_flag){
    
    theta_est <- coefficients(selected_model)['theta']
    if(AV_flag){
      AV_model <- suppressWarnings(suppressMessages(DI(y = y, block = block, density = density, 
                     prop = prop, treat = treat, theta = theta_est, 
                     data = data, DImodel = 'AV')))
      model_list[[length(model_list)+1]] <- AV_model
      names(model_list)[length(model_list)] <- 'AV_model'
    }
    if(FG_flag){
      FG_model <- suppressWarnings(suppressMessages(DI(y = y, block = block, density = density, 
                     prop = prop, treat = treat, theta = theta_est,  
                     data = data, FG = FG, DImodel = 'FG')))
      model_list[[length(model_list)+1]] <- FG_model
      names(model_list)[length(model_list)] <- 'FG_model'
    }
    if(ADD_flag){
      ADD_model <- suppressWarnings(suppressMessages(DI(y = y, block = block, density = density, 
                      prop = prop, treat = treat, theta = theta_est, 
                      data = data, DImodel = 'ADD')))
      model_list[[length(model_list)+1]] <- ADD_model
      names(model_list)[length(model_list)] <- 'ADD_model'
    }
    if(FULL_flag){
      FULL_model <- suppressWarnings(suppressMessages(DI(y = y, block = block, density = density, 
                       prop = prop, treat = treat, theta = theta_est, 
                       data = data, DImodel = 'FULL')))
      model_list[[length(model_list)+1]] <- FULL_model
      names(model_list)[length(model_list)] <- 'FULL_model'
    }
    
    if(!is.na(treat)){
      names(model_list) <- paste0(names(model_list), '_treat')
    }
    names(model_list)[-c(1,2)] <- paste0(names(model_list)[-c(1,2)], '_theta')
    
  } else {
    if(AV_flag){
      AV_model <- suppressWarnings(suppressMessages(DI(y = y, block = block, density = density, 
                     prop = prop, treat = treat, estimate_theta = F,
                     data = data, DImodel = 'AV')))
      model_list[[length(model_list)+1]] <- AV_model
      names(model_list)[length(model_list)] <- 'AV_model'
    }
    if(FG_flag){
      FG_model <- suppressWarnings(suppressMessages(DI(y = y, block = block, density = density, 
                     prop = prop, treat = treat, estimate_theta = F, 
                     data = data, FG = FG, DImodel = 'FG')))
      model_list[[length(model_list)+1]] <- FG_model
      names(model_list)[length(model_list)] <- 'FG_model'
    }
    if(ADD_flag){
      ADD_model <- suppressWarnings(suppressMessages(DI(y = y, block = block, density = density, 
                      prop = prop, treat = treat, estimate_theta = F, 
                      data = data, DImodel = 'ADD')))
      model_list[[length(model_list)+1]] <- ADD_model
      names(model_list)[length(model_list)] <- 'ADD_model'
    } 
    if(FULL_flag){
      FULL_model <- suppressWarnings(suppressMessages(DI(y = y, block = block, density = density, 
                       prop = prop, treat = treat, estimate_theta = F, 
                       data = data, DImodel = 'FULL')))
      model_list[[length(model_list)+1]] <- FULL_model
      names(model_list)[length(model_list)] <- 'FULL_model'
    }
    
    if(!is.na(treat)){
      names(model_list) <- paste0(names(model_list), '_treat')
    }
    
  }
  
  llik <- sapply(model_list, function(x) as.numeric(logLik(x)))
  mAIC <- sapply(model_list, AIC2)
  mAICc <- sapply(model_list, AICc)
  mBIC <- sapply(model_list, BIC2)
  mBICc <- sapply(model_list, BICc)
  ndf <- sapply(model_list, function(x) length(coef(x)))
  
  if (selection == 'Ftest'){
    if (!unique(is.null(FG))){
      message("Since 'Ftest' was specified as selection criterion and functional groups were specified, dropping the ADD model as it is not nested within the FG model.")
      model_list <- model_list[grep("^(?!ADD).*$",names(model_list), value =T, perl = T)]
    }
  }
  selected_model <- switch(selection,
                           Ftest = test_autoDI(model_list = model_list, family = family, treat = treat),
                           AIC = AICsel_autoDI(model_list = model_list, mAIC = mAIC, treat = treat),
                           AICc = AICcsel_autoDI(model_list = model_list, mAICc = mAICc, treat = treat),
                           BIC = BICsel_autoDI(model_list = model_list, mBIC = mBIC, treat = treat),
                           BICc = BICcsel_autoDI(model_list = model_list, mBICc = mBICc, treat = treat)
  )
  
  if(is.null(FG)) {
    message("\nFunctional groups (argument 'FG') were not specified, and therefore not investigated.") 
  }
  if(length(prop) <= 3) {
    message("\nThe 'ADD' variables are only computed for > 3 species cases as the 'ADD' model is not informative for the 2 or 3 species case.") 
  }
  
  message("\nSelected model: ",
          namesub_autoDI(selected_model),
          #        "Formula: ",
          sep = "")
  
  fit_final <- model_list[[selected_model]]
  estimate_theta <- theta_flag
  DImodel <-  fit_final$DIcall$DImodel
  the_final_model <- suppressWarnings(suppressMessages(DI(y = y, prop = prop, block = block,
                                         FG = FG, density = density, treat = treat,
                                         estimate_theta = estimate_theta,
                                         DImodel = DImodel, data = data)))
  the_final_model$DIcall$prop <- prop
  the_final_model$DIcall$DImodel <- DImodel
  the_final_model$DIcall$estimate_theta <- estimate_theta
  the_final_model$DIcall$FG <- FG
  the_final_model$DIcall$treat <- treat
  the_final_model$DIcall$block <- block
  the_final_model$DIcall$density <- density
  the_final_model$DIcall$y <- y
  the_final_model$DIcall$data <- data
  
  return(list('model' = the_final_model, 
              'model_name' = selected_model))
}

autoDI_step3 <- function(selected_model, selection, family) {
  if(length(grep("treat", selected_model$model_name)) == 0) {
    
    message("\n", strrep("-", getOption("width")))
    message("Step 3: No investigation of treatment effect included, since no treatment was specified
        (argument 'treat' omitted)")
    fit_final <- selected_model$model
    selected_model_final <- selected_model$model_name
    
  } else {
    
    fit_treat <- selected_model$model
    fit_notreat <- suppressWarnings(suppressMessages(update_DI(object = fit_treat, treat = NA)))
    
    treat <- fit_treat$DIcall$treat
    
    message("\n", strrep("-", getOption("width")))
    message("Step 3: Investigating the treatment effect\n")
    final_model_list <- list(fit_notreat,
                             fit_treat)
    
    if(length(grep("theta", selected_model$model_name)) == 1) {
      names(final_model_list) <- c(substr(selected_model$model_name, 1, nchar(selected_model$model_name) - 12),
                                   substr(selected_model$model_name, 1, nchar(selected_model$model_name) - 6))
    } else {
      names(final_model_list) <- c(substr(selected_model$model_name, 1, nchar(selected_model$model_name) - 6),
                                   selected_model$model_name)
    }
    
    mAIC_final <- sapply(final_model_list, AIC2)
    mAICc_final <- sapply(final_model_list, AICc)
    mBIC_final <- sapply(final_model_list, BIC2)
    mBICc_final <- sapply(final_model_list, BICc)
    selected_model_final <- switch(selection,
                                   Ftest = test_autoDI(model_list = final_model_list, family = family, treat = treat),
                                   AIC = AICsel_autoDI(model_list = final_model_list, mAIC = mAIC_final, treat = treat),
                                   AICc = AICcsel_autoDI(model_list = final_model_list, mAICc = mAICc_final, treat = treat),
                                   BIC = BICsel_autoDI(model_list = final_model_list, mBIC = mBIC_final, treat = treat),
                                   BICc = BICcsel_autoDI(model_list = final_model_list, mBICc = mBICc_final, treat = treat))
    
    fit_final <- final_model_list[[selected_model_final]]
    
    message("\nSelected model: ",
            namesub_autoDI(selected_model_final), "\n",
            #        "Formula: ",
            sep = "")
  } 
  
  return(list("model" = fit_final,
              "model_name" = selected_model_final))
  
}

autoDI_step4 <- function(prop, data, selected_model, family){
  message("\n", strrep("-", getOption("width")))
  message("Step 4: Comparing the final selected model with the reference (community) model")
  community <- get_community(prop = prop, data = data)
  model_to_compare <- selected_model$model
  model_to_compare_data <- model_to_compare$data
  model_to_compare_data$community <- community
  ref_model <- update(model_to_compare, . ~ . + community,
                      data = model_to_compare_data)
  reftest_autoDI(model_to_compare, ref_model, family)
  
}
