
predict.DI <- function (object, newdata, se.fit = FALSE, type = c("link", 
                                                                  "response"), ...) 
{
  # Make predictions for raw data if no newdata specified
  if (missing(newdata)) {
    if(missing(type)) {
      type <- "link"
    }
    return(predict.glm(object, se.fit = se.fit, type = type, ...))
  }
  
  # Getting species columns and DImodel for adding interactions
  type <- match.arg(type)
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
  prop_cols <- eval(object$DIcall$prop)
  
  if (is.numeric(prop_cols)){
    prop <- names(original_data[, prop_cols])
  } else {
    prop <- prop_cols
  }
  
  
  if (!is.null(prop) & !all(prop %in% colnames(newdata))){
    prop_in_data <- prop[prop %in% colnames(newdata)]
    prop_missing <- prop[!prop %in% prop_in_data]
    if(all(is_near(rowSums(newdata[, prop_in_data]), 1))){
      missing_prop_cols <- matrix(0, ncol = length(prop_missing), nrow = nrow(newdata), dimnames = list(NULL, prop_missing))
      newdata <- cbind(newdata, missing_prop_cols)
      warning(paste0('Species ', paste0(prop_missing, collapse = ','), ' were not present in newdata. These species proportions are assumed to be 0.'))
    } else {
      stop('The species proportions present in newdata don\'t sum to 1.')
    }
  }
  
  if (!is.null(prop) & !all(is_near(rowSums(newdata[, prop]), 1))){
    stop('The species proportions present in newdata don\'t sum to 1.')
  }
  # DI_data doesn't work if dataframe has a single row
  # Adding a dummy row which will be deleted later
  only_one_row <- nrow(newdata) == 1
  if (only_one_row) {
    newdata <- rbind(newdata, rep(0, ncol(newdata)))
  }
  
  # Adding interactions
  theta_flag <- object$coefficients['theta']
  betas <- coef(object)
  if (!DImodel_tag %in% c("ID", "STR")) {
    if (!is.na(theta_flag)){
        theta_value <- coef(object)["theta"]
        betas <- betas[-length(betas)]
    } 
    else {
        theta_value <- 1
    }
    extra_variables <- DI_data(prop = prop, FG = eval(object$DIcall$FG), 
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
          warning(paste0(names(structures[structures == covariate]), ' not supplied in newdata. Calculating for \'', covariate,
                         '\' = ' , 0))
          updated_newdata[, covariate] <- 0
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
          stop(paste0('Values for ', covariate,' given were not present in raw data used for fitting. Predictions can\'t be made for these values.'))
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
    
    extra_data <- model.frame(terms(extra_formula), updated_newdata)
    
    # Matching factors in extra_formula to ones in original_data
    og_factors <- original_data[, sapply(original_data, function(x){is.factor(x) | is.character(x)})]
    common_factors <- intersect(colnames(extra_data), colnames(og_factors))
    
    if (length(common_factors)!=0){
      
      # Levels of all factors in extra_formula
      xlevels <- lapply(common_factors, function(x){levels(original_data[,x])})
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
    
    extra_variables <- DI_data(prop = prop, FG = eval(object$DIcall$FG), 
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
    glm_fit <- glm(fmla, data = original_data_updated, family = object$family)
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
  #print(X_old)
  # glm formula adds NA for non-estimable levles of factors
  # Removing the ones with NA to get only those present in DImodel 
  common <- intersect(names(betas), names(X_old))
  X_new <- as.matrix(X_old[, common])
  
  # Predictions
  y_hat <- as.numeric(X_new %*% betas)

  if (type == "response") {
    inv_link <- object$family$linkinv
    y_hat <- inv_link(y_hat)
  }
  if (se.fit) {
    standard_errors <- as.numeric(sqrt(diag(X_new %*% vcov(object) %*% 
                                              t(X_new))))
    dispersion <- summary(object)$dispersion
    residual.scale <- as.vector(sqrt(dispersion))
    ret_obj <- list(fit = y_hat, se.fit = standard_errors, residual.scale = residual.scale)
  }
  else {
    ret_obj <- y_hat
  }
  return(ret_obj)
}


contrasts_DI <- function(object, contrast, contrast_vars, ...){
  if (missing(object) | class(object)[1]!= 'DI'){
    stop('Please provied a DImodels model object')
  }
  
  if (missing(contrast_vars) & missing(contrast)){
    stop('Provide either one of contrast_vars or constrast')  
  }
  
  if (!missing(contrast_vars) & !missing(contrast)){
    stop('Provide only one of contrast_vars or constrast')  
  }
  
  betas <- coef(object)
  theta_flag <- object$coefficients['theta']
  if (!is.na(theta_flag)){
    theta_value <- coef(object)["theta"]
    betas <- betas[-length(betas)]
  }
  
  if (!missing(contrast_vars)){
    
    if(!is.list(contrast_vars)){
      stop('Contrast variables should be specified as a nested list')
    }
    
    og_data <- object$original_data
    the_C <- matrix(0, ncol = length(betas))
    colnames(the_C) <- names(betas)
    
    for (var in names(contrast_vars)){
      if (!(var %in% colnames(og_data) | var %in% colnames(betas))){
        stop(paste0(var, ' not present in model'))
      }
      
      if(!is.numeric(og_data[,var])){
        if ((length(unlist(contrast_vars[[var]]))/length(unique(og_data[,var]))) %% 1 !=0){
          stop('Lengths of each element of contrasts list should be same as levels of variable in model')
        }
      }
          
      match_names <- paste0(var, levels(og_data[,var]))
      contr.mat <- t(sapply(contrast_vars[[var]], identity))
      colnames(contr.mat) <- match_names
      if (is.null(rownames(contr.mat))){
        rownames(contr.mat) <- paste0(var ,' Test ',1:nrow(contr.mat))
      }
      C_iter <- matrix(0, ncol = length(betas), nrow = nrow(contr.mat))
      colnames(C_iter) <- names(betas)
      common <- intersect(names(betas), match_names)
      C_iter[, common] <- contr.mat[, common]
      rownames(C_iter) <- rownames(contr.mat)
      the_C <- rbind(the_C, C_iter)
    }
    the_C <- the_C[-1,]
    if(class(the_C)[1]=='numeric'){
      the_C <- t(as.matrix(the_C))
      rownames(the_C) <- rownames(contr.mat)
    }
  }
  
  if (!missing(contrast)){
    if(class(contrast)[1] == 'list'){
      if (!all(lengths(contrast) == length(betas))){
        stop('Lengths of each element of contrasts list should be same as number of coefficients in model')
      }
      the_C <- t(sapply(contrast, identity, simplify = TRUE, USE.NAMES = TRUE))
    } else if (class(contrast)[1] == 'numeric'){
      if((length(contrast)/length(betas)) %% 1 !=0){
        stop('Number of elements in contrasts vector should be a multiple of number of coefficients in model')
      }
      the_C <- matrix(contrast, ncol = length(betas), byrow = TRUE)
    } else if (class(contrast)[1] == 'matrix'){
      if(ncol(contrast)!= length(betas)){
        stop('Number of columns in contrast matrix should be same as number of coefficients in model')
      }
      the_C <- contrast
    } else {
      stop('Specify contrast as either a numeric vector, list or matrix')
    }
    if (is.null(rownames(the_C))){
      rownames(the_C) <- paste0('Test ',1:nrow(the_C))
    }
    colnames(the_C) <- names(betas)
  }
  
  cat('Generated contrast matrix:\n')
  print(the_C)
  
  contr.test <- multcomp::glht(object, linfct = the_C, coef = betas, vcov = vcov(object), ...)
  return(contr.test)
}

is_near <- function (x, y, tol = .Machine$double.eps^0.5) {
  abs(x - y) < tol
}