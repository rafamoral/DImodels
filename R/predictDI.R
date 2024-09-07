# Predict function
predict.DI <- function (object, newdata, se.fit = FALSE, 
                        interval = c("none", "confidence", "prediction"), 
                        level = 0.95, weights = 1,
                        type = c("link", "response", "terms"), ...) 
{
  interval <- match.arg(interval)
  
  # Make predictions for raw data if no newdata specified
  if (missing(newdata)) {
    if(missing(type)) {
      type <- "link"
    }
    ret_obj <- predict.glm(object, 
                           se.fit = ifelse(interval != "none", TRUE, se.fit),
                           type = type, ...)
  } else {
    # Getting species columns and DImodel for adding interactions
    type <- match.arg(type)
    newdata <- as.data.frame(newdata)
    DImodel_tag <- object$DIcall$DImodel
    if (is.null(DImodel_tag)){
      DImodel_tag <- 'CUSTOM'
    }
    
    # If custom_formula was used just pass the object to predict.glm to get predictions
    if (DImodel_tag == 'CUSTOM'){
      if(missing(type)) {
        type <- "link"
      }
      return(predict.glm(object, newdata = newdata, se.fit = se.fit, type = type, ...))
    }
    
    original_data <- object$original_data
    model_data <- eval(object$model)
    prop <- attr(object, "prop")
    ID_cols <- attr(object, "ID")
    
    if (!is.null(prop) & !all(prop %in% colnames(newdata))){
      prop_in_data <- prop[prop %in% colnames(newdata)]
      prop_missing <- prop[!prop %in% prop_in_data]
      missing_prop_cols <- matrix(0, ncol = length(prop_missing), nrow = nrow(newdata), dimnames = list(NULL, prop_missing))
      newdata <- cbind(newdata, missing_prop_cols)
      warning(paste0('Species ', paste0(prop_missing, collapse = ', '), ' were not present in newdata. These species proportions are assumed to be 0 and predictions would be returned but might not be meaningful. Check the prop columns in `newdata` to ensure they sum to 1..'))
    } else if (!is.null(prop) & !all(is_near(rowSums(newdata[, prop]), 1))){
      warning('The species proportions present in newdata don\'t sum to 1. \nThe predictions would be returned but might not be meaningful. Check the prop columns in `newdata` to ensure they sum to 1.')
    }
    
    # DI_data doesn't work if dataframe has a single row
    # Adding a dummy row which will be deleted later
    only_one_row <- nrow(newdata) == 1
    if (only_one_row) {
      newdata <- rbind(newdata, newdata)
    }
    
    # Adding interactions
    theta_flag <- object$coefficients['theta']
    betas <- coef(object)
    if (!is.na(theta_flag)){
      theta_value <- coef(object)["theta"]
      betas <- betas[-length(betas)]
    } 
    else {
      theta_value <- 1
    }
    
    if (!DImodel_tag %in% c("ID", "STR")) {
      
      extra_variables <- DI_data(prop = prop, FG = attr(object, "FG"), 
                                 data = newdata, theta = theta_value, what = DImodel_tag)
      if (DImodel_tag == "E") {
        updated_newdata <- data.frame(newdata, E = extra_variables)
      }
      if (DImodel_tag == "AV") {
        updated_newdata <- data.frame(newdata, AV = extra_variables)
      }
      if (DImodel_tag == "ADD") {
        updated_newdata <- data.frame(newdata, extra_variables)
      }
      if (DImodel_tag == "FG") {
        # colnames(extra_variables) <- paste0("FG_", 
        #                                     colnames(extra_variables))
        
        # FG model wasn't working for some reason, so had to assign it this way
        newdata[, 'FG_'] <- extra_variables
        updated_newdata <- newdata
      }
      if (DImodel_tag == "FULL") {
        updated_newdata <- data.frame(newdata, extra_variables, 
                                      check.names = FALSE)
      }
    } else {
      updated_newdata <- newdata
    }
    
    # Grouping ID terms 
    # If not grouping was specified
    if(is.null(ID_cols)){
      ID_cols <- paste0(prop, "_ID")
    }
    
    grouped_IDs <- group_IDs(data = updated_newdata, prop = prop, ID = ID_cols)
    updated_newdata <- cbind(updated_newdata, grouped_IDs)
    
    # Removing the dummy row added
    if (only_one_row) {
      updated_newdata <- updated_newdata[1,]
    }
    
    # Checking for experimental structures
    treat <- eval(object$DIcall$treat)
    density <- eval(object$DIcall$density)
    block <- eval(object$DIcall$block)
    
    structures <- list('treatment' = treat,
                       'density' = density,
                       'block' = block)
    
    for (covariate in structures){
      if (!is.null(covariate)  && !is.na(covariate)){
        # If covariate was supplied as numeric in function call, getting its value
        if (is.numeric(covariate)){
          covariate <- colnames(original_data)[covariate]
        }
        
        if (is.numeric(original_data[, covariate])){
          if ( !(covariate %in% colnames(updated_newdata))){
            warning(paste0(names(structures[structures == covariate]), ' not supplied in newdata. Calculating the prediction for the median value (', median(original_data[, covariate]),') of \'', covariate,
                           '\' from the training data.'))
            updated_newdata[, covariate] <- median(original_data[, covariate])
          }
        } else {
          # Levels of factor covariate in original data
          covariate_levels <- as.factor(unique(original_data[, covariate]))
          # If covariate isn't present in newdata, estimating for base level
          if ( !(covariate %in% colnames(updated_newdata))){
            warning(paste0(names(structures[structures == covariate]), ' not supplied in newdata. Calculating for \'', covariate,
                           '\' = ' , levels(covariate_levels)[1]))
            updated_newdata[, covariate] <- levels(covariate_levels)[1]
          }
          
          # If levels of covariate in newdata not matching ones in original data, stop prediction
          if (! (all(unique(updated_newdata[, covariate]) %in% covariate_levels, na.rm = TRUE))){
            stop(paste0('Values for ', covariate,' given were not present in training data used for fitting. Predictions can\'t be made for these values.'))
          }
          
          # If covariate is supplied as character or numeric, converting to factor
          if (!is.factor(updated_newdata[, covariate])){
            updated_newdata[, covariate] <- factor(updated_newdata[, covariate],
                                                   levels = levels(covariate_levels))
          }
        }
      }
    }
    
    
    # Handling extra formula
    extra_formula <- eval(object$DIcall$extra_formula)
    
    if (! is.null(extra_formula)){
      # If any column from extra_formula is missing in updated_newdata
      e <- try(model.frame(terms(extra_formula), updated_newdata), silent = TRUE)
      if(inherits(e, "try-error")){
        extra_vars <- model.frame(terms(extra_formula), original_data)
        for (covariate in colnames(extra_vars)){
          if(!covariate %in% colnames(updated_newdata)){
            if(is.numeric(extra_vars[, covariate])){
              warning(paste0(names(extra_vars[, covariate]), ' not supplied in newdata. Calculating the prediction for the median value (', median(extra_vars[, covariate]),') of \'', covariate,
                             '\' from the training data.'))
              updated_newdata[, covariate] <- median(extra_vars[, covariate])
            } else {
              # Levels of factor covariate in original data
              covariate_levels <- as.factor(unique(extra_vars[, covariate]))
              # If covariate isn't present in newdata, estimating for base level
              if ( !(covariate %in% colnames(updated_newdata))){
                warning(paste0(names(structures[structures == covariate]), ' not supplied in newdata. Calculating for \'', covariate,
                               '\' = ' , levels(covariate_levels)[1]))
                updated_newdata[, covariate] <- levels(covariate_levels)[1]
              }
              
              # If covariate is supplied as character or numeric, converting to factor
              if (!is.factor(updated_newdata[, covariate])){
                updated_newdata[, covariate] <- factor(updated_newdata[, covariate],
                                                       levels = levels(covariate_levels))
              }
            }
          }
        }
      }
      
      extra_data <- model.frame(terms(extra_formula), updated_newdata)
      
      
      # Matching factors in extra_formula to ones in original_data
      og_factors <- original_data[, sapply(original_data, function(x){is.factor(x) | is.character(x)})]
      common_factors <- intersect(colnames(extra_data), colnames(og_factors))
      
      if (length(common_factors)!=0){
        
        # Levels of all factors in extra_formula
        xlevels <- lapply(common_factors, function(x){levels(as.factor(original_data[,x]))})
        names(xlevels) <- common_factors
        
        for (i in common_factors){
          
          # If levels of factors in extra_formula in newdata not matching ones in original data, stop prediction
          if (! (all(unique(updated_newdata[, i]) %in% xlevels[[i]], na.rm = TRUE))){
            stop(paste0('Values for ', covariate,' given were not present in raw data used for fitting. Predictions can\'t be made for these values.'))
          }
          
          # If factors in extra_formula is supplied as character or numeric, converting to factor
          if (!is.factor(updated_newdata[, i])){
            updated_newdata[, i] <- factor(updated_newdata[,i], levels = xlevels[[i]])
          }
        }
      }
      
      # Having certain functions in extra_formula like poly, ns, bs, etc.
      # cause problems in estimating model.matrix for newdata
      
      # So my solution is to simply refit the model with glm when 
      # such functions are used and then make the prediction
      
      
      fmla <-  eval(object$DIcheck_formula)
      
      extra_variables <- DI_data(prop = prop, FG = attr(object, "FG"), 
                                 data = original_data, theta = theta_value, what = DImodel_tag)
      
      if (DImodel_tag == "E") {
        original_data_updated <- data.frame(original_data, E = extra_variables)
      }
      if (DImodel_tag == "AV") {
        original_data_updated <- data.frame(original_data, AV = extra_variables)
      }
      if (DImodel_tag == "ADD") {
        original_data_updated <- data.frame(original_data, extra_variables)
      }
      if (DImodel_tag == "FG") {
        original_data[, 'FG_'] <- extra_variables
        original_data_updated <- original_data
      }
      if (DImodel_tag == "FULL") {
        original_data_updated <- data.frame(original_data, extra_variables, 
                                            check.names = FALSE)
      } 
      if (DImodel_tag == 'ID' | DImodel_tag == 'STR'){
        original_data_updated <- original_data
      }
      
      # Grouping ID terms 
      # If no grouping was specified
      if(is.null(ID_cols)){
        ID_cols <- paste0(prop, "_ID")
      }
      
      grouped_IDs <- group_IDs(data = original_data_updated, prop = prop, ID = ID_cols)
      original_data_updated <- cbind(original_data_updated, grouped_IDs)
      
      glm_fit <- lm(fmla, data = original_data_updated)
    }
    
    # Calculating response
    # Remove backticks from coefficient names for name-matching
    names(betas) <- gsub(pattern = '`', replacement = '' , x = names(betas))
    
    # glm fmla object from DI_check_and_fit
    if (DImodel_tag == 'CUSTOM'){
      fmla <- eval(object$DIcall$custom_formula)
    } else {
      fmla <-  eval(object$DIcheck_formula)
    }
    
    if (! is.null(extra_formula)){
      Terms <- delete.response(terms(glm_fit))
    } else {
      Terms <- delete.response(terms(formula(fmla)))
    }
    # Model matrix for predictions of newdata
    X_old <- as.data.frame(model.matrix(Terms, data = updated_newdata))
    names(X_old) <- gsub('`','', names(X_old))
    # print(X_old)
    # glm formula adds NA for non-estimable levels of factors
    # Removing the ones with NA to get only those present in DImodel
    common <- intersect(names(betas), names(X_old))
    X_new <- as.data.frame(X_old[, common])
    # 
    # # Predictions
    # y_hat <- as.numeric(X_new %*% betas)
    # 
    # if (type == "response") {
    #   inv_link <- object$family$linkinv
    #   y_hat <- inv_link(y_hat)
    # }
    # if (se.fit) {
    #   standard_errors <- as.numeric(sqrt(diag(X_new %*% vcov(object) %*% 
    #                                             t(X_new))))
    #   dispersion <- summary(object)$dispersion
    #   residual.scale <- as.vector(sqrt(dispersion))
    #   ret_obj <- list(fit = y_hat, se.fit = standard_errors, residual.scale = residual.scale)
    # }
    # else {
    #   ret_obj <- y_hat
    # }
    
    # FG model was failing to give predictions this fixes it
    if(DImodel_tag == "FG" && is.null(extra_formula)){
      if(only_one_row){
        X_new <- rbind(X_new, X_new)
        X_new$FG_ <- extra_variables
        X_new <- X_new[1, ]
      } else {
        X_new$FG_ <- extra_variables  
      }
    }
    
    # Prediction gets messy if we have extra_formula 
    # So manaully making prediction using predict.lm
    if (! is.null(extra_formula)){
      # Terms <- delete.response(terms(glm_fit))
      # predict.lm because this is what is called by predict.glm internally
      ret_obj <- suppressWarnings(predict.lm(glm_fit, newdata = updated_newdata,
                                             se.fit = ifelse(interval != "none", TRUE, se.fit), 
                                             type = ifelse(type == "link", "response", type), ...))
    } else {
      # Terms <- delete.response(terms(formula(fmla)))
      
      ret_obj <- suppressWarnings(predict.glm(object, newdata = if(DImodel_tag == "STR") updated_newdata else X_new, 
                                              se.fit = ifelse(interval != "none", TRUE, se.fit),
                                              type = type, ...))
    }
  }
  
  # Confidence and prediction intervals
  if(object$family$family == "gaussian"){
    if(interval != "none"){
      predictor <- ret_obj$fit
      ip <- ret_obj$se.fit^2
      # Calculating CI/PI
      w <- object$weights
      r <- object$residuals
      rss <- sum(if (is.null(w)) r^2 else r^2 * w)
      df <- object$df.residual
      res.var <- rss/df
      pred.var <- res.var/weights
      
      tfrac <- qt((1 - level)/2, df)
      hwid <- tfrac * switch(interval, confidence = sqrt(ip), 
                             prediction = sqrt(ip + pred.var))
      if (type != "terms") {
        predictor <- cbind(predictor, predictor + hwid %o% 
                             c(1, -1))
        colnames(predictor) <- c("fit", "lwr", "upr")
        ret_obj$fit <- predictor
      }
      # Don't give SE if user didn't ask for it
      if(! se.fit){
        ret_obj$se.fit <- NULL
        ret_obj$residual.scale <- NULL
        ret_obj <- predictor
      }
    }
  }
  
  # Sometimes the prediction could be rank-deficient and lm adds this attribute
  # which could be scary. So drop it
  if(!is.null(attr(ret_obj, "non-estim"))){
    attr(ret_obj, "non-estim") <- NULL
  }
  return(ret_obj)
}

# Contrasts function
contrasts_DI <- function(object, contrast_vars, contrast, ...){
  if (missing(object) | !inherits(object, "DI")){
    stop("Please provied a DImodels model object")
  }
  
  if (missing(contrast_vars) & missing(contrast)){
    stop("Provide either one of `contrast_vars` or `constrast`")
  }
  
  if (!missing(contrast_vars) & !missing(contrast)){
    warning("Provide only one of `contrast_vars` or `constrast`. `contrast_vars` will be ignored.")
    contrast_vars <- NULL
  }
  
  og_data <- object$original_data
  
  # Adjust model coefficients if theta is present
  betas <- coef(object)
  theta_flag <- attr(object, "theta_val")
  if (!is.na(theta_flag)){
    theta_value <- coef(object)["theta"]
    betas <- betas[-length(betas)]
  }
  
  # Branch if contrast_vars are specified
  if (!missing(contrast_vars) && !is.null(contrast_vars)){
    
    # Create contrast matrix according to values specified in contrast_vars
    the_C <- contrast_matrix(object, contrast_vars)
  }
  
  # Branch if contrast is specified
  if (!missing(contrast)){
    if(is.list(contrast)){
      if (!all(lengths(contrast) == length(betas))){
        stop("Lengths of each element of contrasts list should be same as number of coefficients in model")
      }
      the_C <- t(sapply(contrast, identity, simplify = TRUE, USE.NAMES = TRUE))
    } else if (is.matrix(contrast)){
      if(ncol(contrast)!= length(betas)){
        stop("Number of columns in contrast matrix should be same as number of coefficients in model")
      }
      the_C <- contrast
    } else {
      stop("Specify contrast as either a data-frame or matrix.")
    }
    if(is.null(colnames(the_C))){
      colnames(the_C) <- names(betas) 
    }
    the_C <- as.data.frame(the_C)
  }
  
  if (identical(rownames(the_C), as.character(1:nrow(the_C)))){
    rownames(the_C) <- paste0("`Test ", 1:nrow(the_C), "`")
  }
  
  contr_matrix <- as.matrix(the_C)
  
  cat("Generated contrast matrix:\n")
  print(contr_matrix)
  
  contr.test <- multcomp::glht(object, linfct = contr_matrix, coef = betas, vcov = vcov(object), ...)
  return(contr.test)
}

# Helpers for contrast function
add_int_ID <- function(object, newdata){
  # Meta data from model
  prop <- attr(object, "prop")
  ID_cols <- attr(object, "ID")
  DImodel_tag <- attr(object, "DImodel")
  
  # If newdata is not a data-frame convert it to one
  if(!inherits(newdata, "data.frame")){
    newdata <- as.data.frame(newdata)
  }
  
  # DI_data doesn't work if dataframe has a single row
  # Adding a dummy row which will be deleted later
  only_one_row <- nrow(newdata) == 1
  
  if (only_one_row) {
    newdata <- rbind(newdata, newdata)
  }
  
  # Adding interactions
  theta_flag <- attr(object, "theta_val")
  if (!is.na(theta_flag)){
    theta_value <- coef(object)["theta"]
  }
  else {
    theta_value <- 1
  }
  
  if (!DImodel_tag %in% c("ID", "STR")) {
    
    extra_variables <- DI_data(prop = prop, FG = attr(object, "FG"),
                               data = newdata, theta = theta_value, what = DImodel_tag)
    if (DImodel_tag == "E") {
      updated_newdata <- data.frame(newdata, E = extra_variables)
    }
    if (DImodel_tag == "AV") {
      updated_newdata <- data.frame(newdata, AV = extra_variables)
    }
    if (DImodel_tag == "ADD") {
      updated_newdata <- data.frame(newdata, extra_variables)
    }
    if (DImodel_tag == "FG") {
      FG_switch_flags <- c(eval(object$DIcall$treat),
                           eval(object$DIcall$block),
                           eval(object$DIcall$density),
                           eval(object$DIcall$extra_formula))
      # FG model wasn't working for some reason, so had to assign it this way
      if(any(!is.null(FG_switch_flags))){
        colnames(extra_variables) <- paste0("FG_",
                                            colnames(extra_variables))
        updated_newdata <- data.frame(newdata, extra_variables)
      } else {
        newdata[, 'FG_'] <- extra_variables
        updated_newdata <- newdata 
      }
    }
    if (DImodel_tag == "FULL") {
      updated_newdata <- data.frame(newdata, extra_variables,
                                    check.names = FALSE)
    }
  } else {
    updated_newdata <- newdata
  }
  
  # Grouping ID terms
  # If not grouping was specified
  if(is.null(ID_cols)){
    ID_cols <- paste0(prop, "_ID")
  }
  
  grouped_IDs <- group_IDs(data = updated_newdata, prop = prop, ID = ID_cols)
  updated_newdata <- cbind(updated_newdata, grouped_IDs)
  
  # Removing the dummy row added
  if (only_one_row) {
    updated_newdata <- updated_newdata[1,]
  }
  
  return(updated_newdata)
}

contrast_matrix <- function(object, contrast_vars){
  prop <- attr(object, "prop")
  ID <- attr(object, "ID")
  
  # Ensure contrast_vars is specified as a data.frame
  if(!inherits(contrast_vars, "data.frame") && !inherits(contrast_vars, "matrix")){
    warning(paste0("`contrast_vars` should be specified as a <data.frame> or <matrix> containing the contrasts for species proportions in the model, but was specified as a <",
                   class(contrast_vars), ">.\n",
                   "`contrast_vars` will be converted to a <data.frame> but this might not always be possible and might throw errors."))
  }
  contr_data <- as.data.frame(contrast_vars)
  
  # Store additional variables other than proportions or ID effects separately
  # These will be stacked to the final data later
  extra_vars <- names(contr_data)[!names(contr_data) %in% c(prop)]
  
  # The missing proportions in contrast_vars are assumed 0
  prop_missing <- prop[!prop %in% names(contr_data)]
  contr_data[, prop_missing] <- 0
  
  # Data concerning ID effects
  ID_data <- contr_data[, prop]
  # if(!all(is_near(rowSums(ID_data), 0, tol = .Machine$double.eps^0.25))){
  #   warning("The species proportions specified in `contrast_vars` don't all sum to 0 (usually contrasts should sum to 0).\n",
  #           "Assuming this is by choice.")
  # }
  
  # Split species to calculate the net interactions
  positive <- do.call(cbind, sapply(colnames(ID_data), 
                                    function(x) {ifelse(ID_data[, x] >= 0, ID_data[, x], 0)},
                                    simplify = FALSE, USE.NAMES = TRUE))
  negative <- do.call(cbind, sapply(colnames(ID_data), 
                                    function(x) {ifelse(ID_data[, x] < 0, abs(ID_data[, x]), 0)},
                                    simplify = FALSE, USE.NAMES = TRUE))
  
  rownames(positive) <- rownames(negative) <- rownames(as.matrix(ID_data))
  # apply(ID_data, 2, function(x) ifelse(x < 0, abs(x), 0))
  # all.equal(as.data.frame(positive - negative), ID_data)
  
  
  # Calculate interaction terms
  positive <- add_int_ID(object, positive)
  negative <- add_int_ID(object, negative)
  
  # Subtract two components to get back to contrast form
  ID_data <- positive - negative
  
  # Add extra variables onto the matrix
  if(length(extra_vars) > 0){
    # If user has specified any interaction effects trust them and override
    common <- c()
    if(any(extra_vars %in% names(ID_data))){
      # Alert user
      message("Interaction/Identity effects were specified manually in `contrast_vars` ", 
              "using the values specified by user instead of those calculated internally.")
      common <- extra_vars[extra_vars %in% names(ID_data)]
      ID_data[, common] <- contr_data[common]
    }
    the_C <- cbind(ID_data, contr_data[extra_vars[!extra_vars %in% common]])
  } else {
    # If no extra_vars to add then the contrast vector is ready
    the_C <- ID_data
  }
  
  # Special case for FG models
  if(attr(object, "DImodel") == "FG"){
    FG_names <- grep("^FG_.", colnames(the_C), value = TRUE)
    FGs <- the_C[ FG_names]
    colnames(FGs) <- gsub("FG_.", "", colnames(FGs), fixed = TRUE)
    the_C$FG_ <- as.matrix(FGs)
  }
  
  # Add any extra variables in model that weren't specified with value 0
  objTerms <- names(attr(stats::terms(object), "dataClasses"))[-1]
  missing <- objTerms[!objTerms %in% names(the_C)]
  if(length(missing) > 0){
    the_C[, missing] <- 0
  }
  the_C <- as.matrix(model.frame(object$formula[-2], data = the_C))
  return(the_C)
}

# Helpers for predict function
is_near <- function (x, y, tol = .Machine$double.eps^0.5) {
  abs(x - y) < tol
}

# Fortify method for model diagnostics
fortify.DI <- function(model, data = model$model, ...){ 
  # Add proportions to data
  # Only those prop which are different from ID
  prop_idx <- attr(model, "prop")
  IDs <- attr(model, "ID")
  prop_not_in_ID <- prop_idx[!prop_idx %in% IDs]
  prop <- model$data[, prop_not_in_ID]
  data <- cbind(data, prop)
    
  # Add other statistics
  infl <- stats::influence(model, do.coef = FALSE)
  data$.hat <- infl$hat
  data$.sigma <- infl$sigma
  data$.cooksd <- stats::cooks.distance(model, infl)
  data$.fitted <- stats::predict(model)
  data$.resid <- stats::resid(model)
  data$.stdresid <- stats::rstandard(model, infl)
  return(data)
}