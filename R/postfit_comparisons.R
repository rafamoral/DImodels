# Helper function to add missing non-compositional variables to data 
add_extra_vars <- function(object, newdata){
  updated_newdata <- newdata
  original_data <- add_int_ID(object, object$original_data)
  
  # Handling formula
  extra_formula <- eval(object$DIcheck_formula)[-2]
  # browser()
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
          if(is.factor(original_data[, covariate])) 
            covariate_levels <- levels(original_data[, covariate]) 
          else 
            covariate_levels <- levels(as.factor(unique(original_data[, covariate])))
          # If covariate isn't present in newdata, estimating for base level
          if ( !(covariate %in% colnames(updated_newdata))){
            warning(paste0(covariate, ' not supplied in newdata. Calculating for \'', covariate,
                           '\' = ' , covariate_levels[1]))
            updated_newdata[, covariate] <- covariate_levels[1]
          }
          
          # If covariate is supplied as character or numeric, converting to factor
          if (!is.factor(updated_newdata[, covariate])){
            updated_newdata[, covariate] <- factor(updated_newdata[, covariate],
                                                   levels = covariate_levels)
          }
        }
      }
    }
  }
  
  extra_data <- model.frame(terms(extra_formula), updated_newdata)
  
  # browser()
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
        stop(paste0("Values for ", i," given were not present in raw data used for fitting the model.\n",
                    "Predictions can't be made for these values.\nPlease use one of the following ", 
                    paste0(xlevels[i], collapse = ", "), " as values for ", i))
      }
      
      # If factors in extra_formula is supplied as character or numeric, converting to factor
      if (!is.factor(updated_newdata[, i])){
        updated_newdata[, i] <- factor(updated_newdata[,i], levels = xlevels[[i]])
      }
    }
  }
  
  return(updated_newdata)
}

# Helper function to get pairwise combination of every row in a dataframe or 
# against a reference row in the data
get_row_differences <- function(data, ref = NULL) {
  if (!(is.data.frame(data) || is.matrix(data))){
    stop("`data` should be a data.frame or a matrix")
  }
  
  if(! all(sapply(data, is.numeric))){
    stop("All values in `data` should be numeric.")
  }
  
  # browser()
  if (!is.null(ref)) {
    # Validate reference
    if (is.character(ref)) {
      if(! (ref %in% rownames(data))){
        stop("`ref` should be one of the values from the row.names of `data`.")
      }
      ref_row <- data[ref, , drop = FALSE]
      target_rows <- data[setdiff(rownames(data), ref), , drop = FALSE]
    } else if (is.numeric(ref)) {
      if(!(ref >= 1  && ref <= nrow(data))){
        stop(paste0("If specified as a number, `ref` can take values between 1 and `nrow(data)` (", nrow(data),"in this case)."))
      }
      ref_row <- data[ref, , drop = FALSE]
      target_rows <- data[-ref, , drop = FALSE]
    } else {
      stop("ref must be a row name or row index.")
    }
    
    # Differences with just the reference
    diffs <- lapply(rownames(target_rows), function(rn) {
      diff <- as.numeric(data[rn, ]) - as.numeric(ref_row)
      setNames(as.list(diff), colnames(data))
    })
    
    rownames_out <- paste0(rownames(target_rows), "-", rownames(ref_row))
        
  } else {
    # All pairwise combinations
    n <- nrow(data)
    pairs <- combn(n, 2, simplify = FALSE)
    
    diffs <- lapply(pairs, function(idx) {
      diff <- as.numeric(data[idx[1], ]) - as.numeric(data[idx[2], ])
      setNames(as.list(diff), colnames(data))
    })
    
    rownames_out <- sapply(pairs, function(idx) {
      paste0(rownames(data)[idx[1]], "-", rownames(data)[idx[2]])
    })
    
  }
  
  # Combine lists and add proper row-names
  diff_data <- do.call(rbind, lapply(diffs, as.data.frame))
  rownames(diff_data) <- rownames_out
  
  return(diff_data)
}

# Function will accept a DI model object and a data-frame of species
# communities and return all pairwise contrasts along with cld 
compare_communities <- function(object, data, ref = NULL, 
                                adjust = formals(multcomp::adjusted)$type, 
                                verbose = FALSE, alpha.level = 0.05,
                                Letters = c(letters, LETTERS, "."),
                                ...){
  ######################################################################
  # Sanity checking of inputs
  # Ensure object is of class DI
  if (missing(object) | !inherits(object, "DI")){
    stop("Please provide a DImodels model object in `object`.")
  }
  
  # Ensure data is specified as a data.frame
  if (missing(data) | !inherits(data, "data.frame")){
    stop("Please provide a data.frame containing the species communities to be compared in `data`.\n",
         "This data.frame can also contain values for any additional variables in the model.")
  }
  
  # Need more than one commuity for comparison
  if(nrow(data) <= 1){
    stop("The `data` should contain more than 1 community for comparisons.")
  }
  
  # Evaluate if user has not specific any adjustment method
  if (is.name(adjust) || is.call(adjust)) {
    adjust <- eval(adjust)[1]
  }
  
  # Validate specific argument against the available options
  # adjust <- match.arg(adjust)
  
  ############################################################
  # Common variables to be used throughout
  prop <- attr(object, "prop")
  ID <- attr(object, "ID")
  
  # Prepare data for predictions (i.e.,) ensure all necessary columns are present
  data <- as.data.frame(data)
  
  # Any missing species proportions in data are assumed to be 0
  prop_missing <- prop[!prop %in% colnames(data)]
  data[, prop_missing] <- 0
  
  # Ensure species proportions sum to 1
  if(!all(is_near(rowSums(data[, prop]), 1))){
    warning("The species proportions for some communities do not sum to 1.\n",
            "Please check to ensure this is by design.")
  }
  
  # Calculate ID and interaction effects and add them to data
  # Default values for any variables needed for predictions but not 
  # specified in data would also be added here
  # This function creates the entire contrast matrix needed for the cld
  contr_data <- contrast_matrix(object, data)
  
  # If no reference is provided then do all pairwise comparison
  # Get pairwise difference of each row in the contrast matrix
  # All pairwise differences are needed to accurately assign cld
  # But if ref is specified then difference against only that
  # community would be returned
  comp_data <- get_row_differences(contr_data, ref =  ref)
  
  # Get p-values for each pairwise contrast
  comparisons <- contrasts_DI(object = object,
                              contrast = as.matrix(comp_data),
                              verbose = verbose, 
                              ...)

  ## View the tests of comparison
  summary_obj <- summary(comparisons, multcomp::adjusted(type = adjust), ...)               
  summary_obj
  
  # Extract p-values for each test
  test_pvals <- summary_obj$test$pvalues
  names(test_pvals) <- rownames(comp_data)
  
  # This vector of p-values is all we need for assigning the letters.
  # Note that the p-value for each comparison should be named correctly.
  test_pvals  
  
  # If ref is null, letters are assigned to each community
  if(is.null(ref)){
    # Get letters for the communities
    cld_letters <- multcompView::multcompLetters(test_pvals, threshold = alpha.level, 
                                                 Letters = Letters, ...)
    
  # If ref is not null, then we compare against a specific community and 
  # assign letters based on whether that difference is significant or not
  # Note that p-value adjustments for multiple comparisons are accounted 
  # for here as well
  } else {
    cld_letters <- ifelse(test_pvals < alpha.level, Letters[2], Letters[1])
    cld_letters <- c(Letters[1], cld_letters)
    idx <- if(is.character(ref)) which(rownames(data) == ref) else ref
    names(cld_letters) <- c(paste0(rownames(data)[idx], " (ref)"), 
                            rownames(data)[-idx])
  }
  
  return(list("Contrasts" = summary_obj,
              "CLD" = cld_letters))
}


