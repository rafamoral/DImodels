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
            stop(paste0('Values for ', i,' given were not present in raw data used for fitting. Predictions can\'t be made for these values.'))
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
  prop <- model$data[, prop_not_in_ID, drop = FALSE]
  data <- cbind(data, prop)
  
  # Add factor version of any categorical variables in the model
  cols_to_add <- setdiff(colnames(model$original_data), 
                         colnames(data))
  data <- cbind(data, model$original_data[, cols_to_add])
  
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