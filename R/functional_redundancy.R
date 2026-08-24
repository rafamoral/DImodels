# Helper function to create a class of delta_sd. Will be used to calculate 
# equivalence margin as a multiplier of model residual sd
delta_sd <- function(k = 0.5) {
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k <= 0) {
    stop("Please specify `k` as a single positive number.")
  }
  
  structure(k, class = "delta_sd")
}

# Helper function to get residual model sd. Very simple now but may be extended in future
get_model_sigma <- function(obj){
  sig_val <- sigma(obj)
  return(sig_val)
}

# Helper function for pretty formatting of p-values
format_p <- function(p) {
  if (p < 0.0001) {
    "< 0.0001"
  } else if (p < 0.01) {
    paste0("= ", formatC(p, format = "f", digits = 4))
  } else if (p < 0.1) {
    paste0("= ", formatC(p, format = "f", digits = 3))
  } else {
    paste0("= ", formatC(p, format = "f", digits = 2))
  }
}

# Helper functions for situations when someone uses dplyr operations to modify
# output of fr_communities
# Code adapted from example at https://www.r-bloggers.com/2023/04/extending-data-frames/
fr_communities_reconstruct <- function(x, to) {
  attributes(x)$redundant <- attr(to, "redundant")
  attributes(x)$ID_cols <- attr(to, "ID_cols")
  attributes(x)$class <- attr(to, "class")
  
  x
}

`[.fr_communities` <- function(x, ...) {
  
  out <- NextMethod("[")
  
  if (is.data.frame(out)) {
    out <- fr_communities_reconstruct(out, x)
  }
  out
}

`names<-.fr_communities` <- function(x, value) {
  out <- NextMethod("names<-")

  fr_communities_reconstruct(out, x)
}

dplyr_reconstruct.fr_communities <- function(data, template) {
  fr_communities_reconstruct(data, template)
}

# Registration for dplyr_reconstruct S3 method. Will register only if
# dplyr is loaded in namespace
.onLoad <- function(libname, pkgname) {
  register <- function(...) {
    if (isNamespaceLoaded("dplyr")) {
      registerS3method(
        "dplyr_reconstruct", "fr_communities",
        dplyr_reconstruct.fr_communities,
        envir = asNamespace("dplyr")
      )
    }
  }

  # If dplyr already loaded when DImodels loads
  register()

  # If dplyr loaded later in the same session
  setHook(packageEvent("dplyr", "onLoad"), register)
  invisible()
}

# Function to extract relevant set of communties to compare in a objectl
fr_communities <- function(prop, redundant){
  
  # If prop is empty, try to retrieve it from fitted objectl.
  # Throw error if not possible
  if(missing(prop)){
    stop(paste("`prop` should be a vector specifying column names identifying the species proportions."))
  }
  
  if(missing(redundant)){
    stop(paste("`redundant` should be a character vector containing species names that need to be tested for functional redundancy."))
  }
  
  if(length(unique(redundant)) < 2){
    stop("There should be at least two unique species specified in `redundant`")
  }
  
  if(any(!redundant %in% prop)){
    stop(paste0("All names specified in `redundant` should be present in `prop`. ",
               paste0(redundant[!redundant %in% prop], collapse = ", "),
               " is not present in `prop`.\n",
               "`prop` had values ", paste0(prop, collapse = ", ")))
  }
  
  # Set of redundant species
  redundant <- unique(redundant) # Incase user repeats same species
  R_set <- redundant
  # Complete of R_set
  Rp_set <- prop[!prop %in% redundant]
  
  # To compare monocultures of species in R_set
  monos <- as.data.frame(diag(1, length(R_set)))
  names(monos) <- R_set
  
  # Compare monocultures in R_set with their equi-proportional community
  within_comms <- rbind(monos, rep(1/length(R_set), times = length(R_set)))
  within_comms <- cbind(within_comms, 
                        matrix(0, 
                               nrow(within_comms), 
                               ncol = length(Rp_set),
                               dimnames = list(NULL, Rp_set)))
  within_comms$Set <- c(rep("Monos", times = length(R_set)), "Within")
  within_comms$Comm_ID <- c(R_set, "equi")
  within_comms <- within_comms[, c("Set", "Comm_ID", prop)]
  rownames(within_comms) <- c(R_set, "equi")
  
  # Binary mixtures of each species in Rp_set with all in R_set
  pairs <- expand.grid(R_set, Rp_set)

  between_comms <- as.data.frame(
    outer(pairs$Var1, prop, `==`) +
    outer(pairs$Var2, prop, `==`)
  )
  between_comms <- between_comms/rowSums(between_comms)
  colnames(between_comms) <- prop
  between_comms$Set <- rep("Between", times = nrow(between_comms))
  between_comms$Comm_ID <- pairs$Var2
  
  # Add rownames only if there are between communities
  if(nrow(between_comms) > 0){
    rownames(between_comms) <- paste0(pairs$Var2, "_", pairs$Var1)
  } 
  
  # Merge the data.frames
  ret_df <- rbind(within_comms, between_comms)
  
  # Adding some attributes to the data. Helpful for operations later
  attr(ret_df, "redundant") <- redundant
  
  attr(ret_df, "ID_cols") <- prop
  
  class(ret_df) <- c("fr_communities", class(ret_df))
  
  return(ret_df)
}

# This function accepts output from fr_communities and creates the 
# contrast matrix with the specific pairwise contrasts needed to confirm
# functional redundancy
fr_contrasts <- function(data, object, 
                         int_cols){
  if(missing(object)){
    stop("Please specify a regression model in `object`")
  }
  
  if(missing(data)){
    stop("`data` is missing.\n",
         "Specify the data.frame returned by the `fr_communities()` function ", 
         "after adding in all relevant terms needed for model predictions.")
  }
  
  # Species which are tested to be redundant
  redundant <- attr(data, "redundant")
  if(is.null(redundant) || !inherits(data, "fr_communities")){
    stop("`data` should be a data.frame returned by the `fr_communities()` function. ",
         "\nIf you get this error even after using `fr_communities()`, please review ",
         "your data manipulation pipeline as it may have corrupted some attributes.")
  }
  
  if(is.null(data$Set)){
    stop("The `Set` column in the data is important. Please rerun the function with it included.")
  }
  
  if(is.null(data$Comm_ID)){
    stop("The `Comm_ID` column in the data is important. Please rerun the function with it included.")
  }
  
  # If int_cols is missing try to salvage them from data
  if(missing(int_cols) && !inherits(object, "DI")){
    stop("`int_cols` is mandatory. Specify the names or indices of the columns in `data` containing interactions as a vector.",
         "If there are no interaction columns in model, specify `int_cols` as NULL.")
  }
  
  # Sanity checks for int_cols
  if(!missing(int_cols) && !is.null(int_cols)){
    if (is.character(int_cols) && any(!int_cols %in% names(data)))
      stop("Can't select columns: ", 
           paste(int_cols[!int_cols %in% names(data)], collapse = ", "),
           "\n as they do not exist in `data`.",
           call. = FALSE)
    if (is.numeric(int_cols) && any(!int_cols %in% seq_len(ncol(data))))
      stop("Can't select columns past the end: ", 
           "\nThere are only ", ncol(data), " columns, but ",
           "positions beyond ", ncol(data), " were specified.",
           call. = FALSE)
  } 
  
  
  # If object is DImodel, add interactions if they are not specified by user. 
  if(inherits(object, "DI")){
    # This will only add IDs and interactions if user didn't specify them initially.
    # If user specified values are present they are left unchanged.
    data <- add_int_ID(object, data)
    
    # If int_cols is missing try to salvage them from data
    if(missing(int_cols)){
      int_cols <- attr(data, "int_cols")
    }
  }

  # Position of within index needed for adding extra row with 0 interactions
  within_idx <- which(data$Set == "Within")
  int_0_row <- data[within_idx, , drop = FALSE]
  rownames(int_0_row) <- "mono_av"
  int_0_row$Comm_ID <- "mono_av"
  # For mono average data make interaction 0 (if they're present)
  if(!is.null(int_cols)){
    int_0_row[, int_cols] <- 0
  }
  
  # For the within test row create copy with interactions as zero
  data <- rbind(data[seq_len(within_idx), ],
                       int_0_row,
                       data[-seq_len(within_idx), ])
  
  
  # First try to create model matrix normally. If this works fine, move straight ahead
  mm_err <- NULL
  mod_mat <- tryCatch(model.matrix(formula(object)[-2], data),
                      error = function(e) { mm_err <<- e; NULL })
  # If model.matrix fails, and model is fit using DImodels fallback to contrast_matrix()
  if(is.null(mod_mat) && inherits(object, "DI")){
    # Try constructing contrast matrix using DI setup. 
    mod_mat <- tryCatch(suppressMessages(contrast_matrix(object, data)), 
                        error = function(e) NULL)
  }
  
  # If model.matrix failed and model wasn't fit with DImodels or if even 
  # contrast_matrix fails. Throw original error
  if(is.null(mod_mat)){
    # Check that all variables required by the model formula are present
    required <- all.vars(formula(object)[-2])
    missing <- setdiff(required, names(data))
    
    # Custom error message if any variable is missing
    if (length(missing) > 0) {
      stop(
        "The following variables required by the model formula ",
        "are missing from `data`: ",
        paste("\n  -", missing, collapse = ""),
        call. = FALSE
      )
    }
    
    # Checking if any factor variable in model was specified as a character string
    expected_factors <- intersect(names(object$xlevels), names(data))
    wrong_type <- expected_factors[
      !sapply(data[expected_factors], is.factor)
    ]
    
    if (length(wrong_type) > 0) {
      stop(
        "The following variables were fitted as factors in the model but are not factors in `data`: ",
        paste("\n  -", wrong_type, collapse = ""),
        "\nConvert them to factors having same levels as in raw data using `factor(values, levels = c(...))`.",
        call. = FALSE
      )
    }
    
    # Fallback generic error message from model.matrix
    stop("Failed to construct model matrix from `data`: ",
         conditionMessage(mm_err), call. = FALSE)
  }
  
  rownames(mod_mat) <- rownames(data)
  
  monos_idx <- which(data$Set == "Monos")
  within_idx <- which(data$Set == "Within")
  between_idx <- which(data$Set == "Between")
  
  equiv_groups <- ifelse(data$Set == "Between", data$Comm_ID, data$Set)
  
  contr_list <- lapply(unique(equiv_groups), function(x){
    ret_mat <- get_row_differences(mod_mat[equiv_groups == x, , drop = FALSE])
    ret_mat
  })
  
  idxs <- sapply(contr_list, nrow)
  names(idxs) <- unique(equiv_groups)
  contr_mat <- do.call(rbind, contr_list)
  
  set_idxs <- c(idxs[c("Monos", "Within")],
                "Between" = sum(idxs[!names(idxs) %in% c("Monos", "Within")]))
  attr(contr_mat, ".Equiv_group") <- rep(names(idxs), idxs)
  attr(contr_mat, "Set") <- rep(names(set_idxs), set_idxs)
  attr(contr_mat, "redundant") <- redundant
  class(contr_mat) <- c("fr_contrasts", class(contr_mat))
  
  return(contr_mat)
}

# This function accept the fr_contrast matrix and performs the equivalence tests
fr_tost <- function(data, object,  
                    delta = delta_sd(0.5), alpha.level = 0.05) {
  

  if(missing(object)){
    stop("Please specify a regression model in `object`")
  }
  
  if(missing(data)){
    stop("`data` is missing.\n",
         "Specify the data.frame containing pairwise contrasts returned by the `fr_contrasts()` function.")
  }
  
  if(!inherits(data, "fr_contrasts")){
    stop("`data` should be a data.frame returned by the `fr_contrasts()` function.")
  }
  
  if(!inherits(data, "matrix")){
    if(inherits(data, "data.frame")){
      cont_mat <- as.matrix(data)
      attr(cont_mat, "Set") <- attr(data, "Set")
      attr(cont_mat, ".Equiv_group") <- attr(data, ".Equiv_group")
      attr(cont_mat, "redundant") <- attr(data, "redundant")
    } else {
      stop("`data` should be a matrix.")
    }
  }
  
  # Resolving delta
  ## If delta is missing, warning should be thrown for user
  if (missing(delta)) {
    warning(
      "No `delta` supplied; using '0.5 x model residual SD' as the default equivalence margin.",
      "\nWe suggest specifying `delta` explicitly with values based on domain knowledge for better inference.")
  }
  
  # Delta will either be specified as delta_sd(k)
  if (inherits(delta, "delta_sd")) {
    delta_str <- paste0("delta_sd(", round(unclass(delta), 2), ")") # Used later in printing
    delta <- unclass(delta) * get_model_sigma(object)
  } else {
    # Or as a numeric vector of length 1 or 2, in which case it needs to be validated
    if (!is.numeric(delta) ||
        length(delta) > 2 || length(delta) < 1 || 
        !all(is.finite(delta))) {
      stop(
        "Please specify `delta` as a positive finite vector of length one or two ",
        "defining equivalence margin based on domain knowledge or specify it using", 
        "`delta_sd(k)` where k is a positive number.")
    }
    
    if(length(delta) == 1 && any(delta < 0)){
      stop(
        "If specifying a single value for `delta`, it should be a positive number."
      )
    }
    
    if(length(delta) == 2 && (all(delta < 0) || all(delta > 0))){
      stop(
        "If specifying `delta` as a two-element vector, the first value should be negative while second should be positive."
      )
    }
    delta_str <- "User specified range"
  }
  
  # If single value specified for delta then both lower and higher bounds are same  (i.e., symmetric interval)
  if(length(delta) == 1){
    delta_lo = -abs(delta)
    delta_hi = abs(delta)    
  } else {
    # Specifying delta as c(low, hi) allows to have different bounds (i.e., asymmetric interval)
    delta_lo = delta[1]
    delta_hi = delta[2]
  }
  
  if (any(delta_hi <= delta_lo)) {
    stop("`delta_lo` must be lower than `delta_hi`.")
  }
  
  # Model objects
  betas <- coef(object)
  V_mat <- vcov(object)
  # If DI model was fitted with theta, it needs to be removed from coefficients 
  # for the computations here
  if(inherits(object, "DI") && !is.na(attr(object, "theta_val"))){
    betas <- betas[names(betas) != "theta"]
  }
  
  if (ncol(cont_mat) != length(betas)) {
    stop(paste("The number of columns in `data` should be same as number of coefficients in model.",
               "`data` has", ncol(cont_mat),"columns but object has", length(betas), "coefficients."))
  }
  
  # Contrast estimate
  nas <- is.na(betas)
  
  if (any(cont_mat[, nas, drop = FALSE] != 0)) {
    stop("One or more contrasts depend on a non-estimable model coefficient.",
         "Check the fitted model for overparmeterisation.", call. = FALSE)
  }
  
  C <- cont_mat[, !nas, drop = FALSE]
  b <- betas[!nas]
  V <- V_mat[!nas, !nas, drop = FALSE]
  
  est <- as.numeric(C %*% b)
  if (anyNA(est)) {
    warning("One or more contrasts produced estimates of NA.",
            "Compare `data` columns against model coefficients.")
  }
  
  # Faster way to compute diag(CVC')
  var_est <- rowSums((C %*% V) * C) 
  # Contrast SE
  se <- sqrt(var_est)
  df <- df.residual(object)

  # Regular t-test comparing if contrast significantly differs from zero
  # Not relevant for equivalence, computed only as an additional artifact
  t_raw <- est / se
  p_raw <- 2 * pt(-abs(t_raw), df)
  
  # TOST
  t_lo <- (est - delta_lo) / se     # H0: true <= delta_lo
  t_hi <- (est - delta_hi) / se     # H0: true >= delta_hi
  p_lo <- pt(t_lo, df, lower.tail = FALSE)
  p_hi <- pt(t_hi, df, lower.tail = TRUE)
  # TOST is intersection union test, thus Type-I error capped at max of p1 and p2.
  p_tost <- pmax(p_lo, p_hi)
  
  crit <- qt(1 - alpha.level, df)
  ci_level <- 1 - 2 * alpha.level
  ci_lo <- est - crit * se
  ci_hi <- est + crit * se
  
  equivalent <- p_tost < alpha.level
  
  # Community identifiers for each contrast
  comm_names <- strsplit(rownames(cont_mat), "-")
  comm1 <- sapply(comm_names, `[`, 1)
  comm2 <- sapply(comm_names, `[`, 2)

  ret_df <- data.frame(
    Set         = attr(cont_mat, "Set"),
    Equiv_group = attr(cont_mat, ".Equiv_group"),
    Comm1        = comm1,
    Comm2        = comm2,
    Estimate     = est,
    SE           = se,
    df           = df,
    alpha        = alpha.level,
    CI_level     = ci_level,
    CI_lo        = ci_lo,
    CI_hi        = ci_hi,
    delta_lo     = delta_lo,
    delta_hi     = delta_hi,
    t_raw        = t_raw,
    p_raw        = p_raw,
    p_lower      = p_lo,
    p_upper      = p_hi,
    p_tost       = p_tost,
    Equivalent   = equivalent
  )
  
  # Assign specific class to help with pretty printing later
  class(ret_df) <- c("fr_tost", class(ret_df))
  attr(ret_df, "delta_str") <- delta_str
  attr(ret_df, "redundant") <- attr(cont_mat, "redundant")
  attr(ret_df, "contrast_matrix") <- cont_mat
  
  return(ret_df)
}


# Helper function to perform functional redundancy test given
# output of fr_communities with interaction columns
fr_test <- function(data, object, int_cols, 
                    delta = delta_sd(0.5), alpha.level = 0.05){
  if(missing(object)){
    stop("Please specify a regression model in `object`")
  }
  
  if(missing(data)){
    stop("`data` is missing.\n",
         "Specify the data.frame returned by the `fr_communities()` function ", 
         "after adding in all relevant terms needed for model predictions.")
  }
  
  ## If delta is missing, warning should be thrown for user
  if (missing(delta)) {
    warning(
      "No `delta` supplied; using '0.5 x model residual SD' as the default equivalence margin.",
      "\nWe suggest specifying `delta` explicitly with values based on domain knowledge for better inference.")
  }
  
  # If int_cols is missing try to salvage them from data
  if(missing(int_cols) && !inherits(object, "DI")){
    stop("`int_cols` is mandatory. Specify the names or indices of the columns in `data` containing interactions as a vector.",
         "If there are no interaction columns in model, specify `int_cols` as NULL.")
  }
  
  # Generate contrast matrix
  contr_mat <- fr_contrasts(data = data, 
                            object = object,
                            int_cols = int_cols) 
  
  # Two one side t-tests results for contrast matrix
  tost_results <- fr_tost(data = contr_mat,
                          object = object, 
                          delta = delta,
                          alpha.level = alpha.level) 

  return(tost_results)
}

# functional redundancy wrapper for DImodels package.
# Only a fitted DImodels and redundant vector are needed.
fr_test_DI <- function(object, redundant, 
                       delta = delta_sd(0.5),
                       alpha.level = 0.05){
  
  # objectl can't be missing
  if(missing(object)){
    stop("Specify a fitted DI model object in `object`.")
  }

  # Flag is objectl is fit using DIobjectls
  is_DI <- inherits(object, "DI")
  if(isFALSE(is_DI)){
    stop("Please provide a model fitted using the `DImodels` package in `object`.")
  }
  
  if(missing(redundant)){
    stop(paste("`redundant` should be a character vector containing species names that need to be tested for functional redundancy."))
  }
  
  ## If delta is missing, warning should be thrown for user
  if (missing(delta)) {
    warning(
      "No `delta` supplied; using '0.5 x model residual SD' as the default equivalence margin.",
      "\nWe suggest specifying `delta` explicitly with values based on domain knowledge for better inference.")
  }
  
  # Species proportions
  prop <- attr(object, "prop")
  ID_vals <- attr(object, "ID")
    
  communities <- fr_communities(prop = prop, 
                                redundant = redundant)
  
  # Generate contrast matrix
  contr_mat <- suppressWarnings(
    fr_contrasts(data = communities, 
                 object = object)
    ) 
  
  # Two one side t-tests results for contrast matrix
  tost_results <- fr_tost(data = contr_mat,
                          object = object, 
                          delta = delta,
                          alpha.level = alpha.level) 

  return(tost_results)
}

# Helper function to print the results of functional redundancy test nicely
print.fr_tost <- function(x, ...){
  
  if(!inherits(x, "fr_tost")){
    stop("`x` should be a data.frame returned by the `fr_tost()` function.")
  }
  
  data <- x
  redundant <- attr(data, "redundant")
  redundancy_str <- if (length(redundant) == 2) {
    paste(redundant, collapse = " and ")
  } else {
    paste0(
      paste(redundant[-length(redundant)], collapse = ", "),
      ", and ",
      redundant[length(redundant)]
    )
  }
  
  alpha <-  unique(data$alpha)
  delta_lo <-  unique(data$delta_lo)
  delta_hi <-  unique(data$delta_hi)
  
  if(attr(data, "delta_str") == "User specified range") 
    delta_text <- "as defined by user" 
  else 
    delta_text <- paste0("quantified using ", 
                         gsub("[^0-9.]", "", attr(data, "delta_str")),
                         " times residual SD of model") 
  
  
  monos_data <- data[data$Set == "Monos", ]
  within_data <- data[data$Set == "Within", ]
  between_data <- data[data$Set == "Between", ]
  
  overall_p <-  max(data$p_tost)
  monos_p <- max(monos_data$p_tost)
  within_p <- max(within_data$p_tost)
  # If all species are tested for redundancy then there's be no between contrasts
  between_p <- if(nrow(between_data) > 0) max(between_data$p_tost) else NA
  
  con_width <- getOption("width", 80)
  cat(strrep("-", con_width), "\n")
  cat(sprintf("\nTesting functional redundancy for species %s.\n\n", redundancy_str))
  cat(strrep("-", con_width), "\n")
  
  cat(sprintf("\nEquivalence tests conducted at the alpha = %s level (i.e, using %s CIs) with a delta margin of (%s, %s) %s.\n\n",
              alpha, paste0(100*(1-2*alpha), "%"), 
              format(round(delta_lo, 2), trim = TRUE),
              format(round(delta_hi, 2), trim = TRUE),
              delta_text))
  
  cat(sprintf("Null hypothesis (H0): The species are not functionally redundant.\n"))
  cat(sprintf("Alternative hypothesis (H1): The species are functionally redundant.\n"))
  cat(paste0("\n", strrep("-", con_width), "\n"))
  cat(sprintf("\nResults: \n"))
  
  if(isTRUE(all(data$Equivalent))){
    
    cat(sprintf("\nThe species %s are functionally redundant (p %s).\n\n",
                redundancy_str, format_p(overall_p)))
    
    success_text <- paste0("All equivalence criteria were satisfied.\n\n",
                           sprintf("  - Monocultures were found to have equivalent performance (p %s). \n", format_p(monos_p)),
                           sprintf("  - Within-group interactions were found to be equivalent to 0 (p %s). \n", format_p(within_p)),
                           sprintf("  - Between-group interactions with all other species were found to be equivalent (p %s).\n", format_p(between_p)))
    cat(success_text)
    cat(paste0("\n", strrep("-", con_width), "\n"))
  } else {
    cat(sprintf("\nFunctional redundancy was not established for %s (p %s) due to the following reasons.\n\n",
                redundancy_str, format_p(overall_p)))
    
    fail_reasons <- ""
    
    if(!isTRUE(all(monos_data$Equivalent))){
      fail_reasons <- paste0(fail_reasons, 
                             "  - Monoculture performances",
                             " were not found to be equivalent (p ", format_p(monos_p), ").\n")
    }
    
    if(!isTRUE(all(within_data$Equivalent))){
      fail_reasons <- paste0(fail_reasons, 
                             "  - Within-group interactions", 
                             " were not found to be equivalent to 0 (p ", format_p(within_p), ").\n")
    }
    
    if(!isTRUE(all(between_data$Equivalent))){
      fail_reasons <- paste0(fail_reasons, 
                             "  - Between-group interactions", 
                             " were not found to be equivalent for the following (p ", 
                             format_p(between_p), "):\n")
      sets <- unique(between_data$Equiv_group)
      for (s in sets){
        # If any pairwise contrast is not equivalent, print it
        if(any(!between_data[between_data$Equiv_group == s, "Equivalent"])){
          p_val <- max(between_data[between_data$Equiv_group == s, "p_tost"])
          fail_reasons <- paste0(fail_reasons,
                                 "    * ", sub(".*with ", "", s), 
                                 " (p ", format_p(p_val), ")\n")
        }
      }
    }
    cat(fail_reasons)
    cat(paste0("\n", strrep("-", con_width), "\n"))
    
  }
  
  cat("Note: All tests follow the intersection-union principle and the respective overall p-value is the maximum p-value across all lower level constituent tests.\n")
  
  return(invisible(data))
}

# Visualise results of functional redundancy test
plot.fr_tost <- function(x,
                         xlab = "Difference in predictions",
                         ...) {
  
  if(!inherits(x, "fr_tost")){
    stop("`x` should be a data.frame returned by the `fr_tost()` function.")
  }
  
  data <- x
  # Equivalence margin boundaries
  delta_lo <- unique(data$delta_lo)
  delta_hi <- unique(data$delta_hi)
  
  # Set should be present in the data for grouping the contrasts
  # If absent, add a generic Set = "All"
  if (!"Set" %in% names(data)) {
    data$Set <- "All"
  }
  
  # Contrast groups in the data
  # Assuming user didn't change anything, 
  # should be only monos, within, and between
  groups <- unique(data$Set)
  
  criterion_labels <- c(
    Monos   = "Monoculture differences",
    Within  = "Within-group interactions",
    Between = "Between-group interactions"
  )
  
  # If groups match desired set, give them informative labels,
  # else leave them unchanged
  group_labels <- ifelse(
    groups %in% names(criterion_labels),
    criterion_labels[groups],
    groups
  )
  names(group_labels) <- groups
  
  # Y position of each contrast. Leave 1 blank space for each separator (group label)
  # Positions are reversed so labels start from top to bottom
  # Numbering starts from # rows in data + # groups and counts down to 1
  group_positions <- sapply(seq_along(groups), 
                            function(x) {
                              1 + nrow(data) + length(groups) - (which(data$Set == groups[x]) + x)
                              })
  names(group_positions) <- groups
  
  # Add new positions to data
  data$y <- NA_integer_
  for (g in names(group_positions)) {
    rows <- which(data$Set == g)
    data$y[rows] <- group_positions[[g]]
  }
  
  # Y positions for the separator (the blank spaces left in group_positions)
  separator_y <- sapply(groups, 
                        function(x) max(group_positions[[x]]) + 1)
  
  # Plot x-axis range
  x_range <- range(
    c(
      data$CI_lo,
      data$CI_hi,
      delta_lo,
      delta_hi
    ),
    finite = TRUE
  )
  
  # 15% padding on padding on left side of plot
  x_pad <- diff(x_range) * 0.15
  
  # Failsafe is x_range is zero
  if (x_pad == 0) {
    x_pad <- 1
  }
  
  xlim <- c(
    x_range[1] - x_pad,
    x_range[2]
  )
  
  # y is always continuous, i.e., 1:nrow(data) + length(groups),
  # padding is 0.75 data units on both ends
  ylim <- c(
    min(data$y) - 0.75,
    max(data$y) + 0.75
  )
  
  contrast_cex <- 0.8 # Font size for the contrast labels
  heading_cex  <- 0.9 # Font size for the group heading labels
  
  # Adjusting left margin on plot to allow labels to display nicely
  # Values chosen with trial and error
  row_label_w <- max(strwidth(rownames(data), units = "inches", cex = contrast_cex))
  heading_w   <- max(strwidth(paste(group_labels, "(p 0.000)"),
                              units = "inches", cex = heading_cex))
  
  line_in <- par("csi")
  heading_line <- 1 + row_label_w / line_in + 0.4   # +0.4 line gap between tiers
  
  left_margin <- heading_line * line_in + heading_w - 0.5
  
  # Edit plot parameters but restore to user parameters on exit
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  
  # Plot margins
  par(
    mai = c(
      1.5,             # bottom
      left_margin,     # left
      0.8,             # top
      0.5              # right
    )
  )
  
  # Plot blank canvas
  plot(
    NA,
    xlim = xlim, ylim = ylim,
    yaxt = "n",  # Suppress plotting y-axis
    xlab = xlab, ylab = "",
    bty = "n"    # Suppress box around plot boundaries
  )

  # Equivalence region
  rect(
    xleft = delta_lo,
    ybottom = ylim[1],
    xright = delta_hi,
    ytop = ylim[2],
    col = adjustcolor(
      "steelblue1",
      alpha.f = 0.6
    ),
    border = NA
  )
  
  # Equivalence boundaries
  abline(
    v = c(delta_lo, delta_hi),
    lty = 2,
    lwd = 1.25,
    col = "steelblue4"
  )
  
  # Null line
  abline(
    v = 0,
    lty = 1,
    lwd = 1.5
  )
  
  # Confidence intervals
  segments(
    x0 = data$CI_lo,
    y0 = data$y,
    x1 = data$CI_hi,
    y1 = data$y,
    lwd = 2.5
  )
  
  # Point estimates
  points(
    data$Estimate,
    data$y,
    pch = 19,
    cex = 1.35
  )
  
  # Row labels
  axis(
    side = 2,
    at = data$y,
    labels = rownames(data),
    las = 1,
    cex.axis = 0.8,
    tick = FALSE
  )
  
  # Equivalence margin annotation
  mtext(
    sprintf(
      "Equivalence margin: (%s, %s)",
      format(delta_lo, digits = 2),
      format(delta_hi, digits = 2)
    ),
    side = 3,   # On top of plot
    line = 0.1, # On margin line 1 
    at = 0,     # At x = 0
    cex = 0.9,
    col = "steelblue4"
  )
  
  # p-values for individual contrasts
  text(
    data$Estimate,
    data$y + 0.25,
    labels = paste0(
      "p ",
      sapply(data$p_tost, format_p)
    ),
    cex = 0.8,
    pos = 4
  )
  
  # Horizontal separators between community sets
  abline(
    h = separator_y,
    col = "grey60",
    lwd = 1
  )
  
  # Labels for each set
  mtext(
    sapply(names(group_labels), function(x){
      idx <- which(data$Set == x)
      sprintf(
        "%s (p %s)",
        group_labels[x],
        format_p(max(data[idx, "p_tost"], na.rm = T))
      )
    }),
    side = 2,
    at = separator_y,
    line = 0.5,
    las = 1,
    font = 2,
    cex = 0.9,
    col = "grey30"
  )
  
  # Overall p-value for functional redundancy test
  overall_p <- max(
    data$p_tost,
    na.rm = TRUE
  )
  mtext(
    sprintf(
      "Overall test p-value %s",
      format_p(overall_p)
    ),
    side = 3,  # add text on top of plot
    line = 1,  # On margin line 1 
    at = 0,    # at x = 0
    cex = 0.9
  )
  
  # Footnote
  mtext(
    paste(
      "Note: All tests follow the intersection-union principple and thus,",
      "overall p-values are the maximum p-value across lower level",
      "constituent equivalence tests."
    ),
    side = 1,
    line = 4,
    cex = 0.65,
    col = "grey40"
  )
  
  invisible(NULL)
}

# Summary function for additional information about test
summary.fr_tost <- function(object,
                            max.print = 20,
                            ...) {
  
  if(!inherits(object, "fr_tost")){
    stop("`object` should be a data.frame returned by the `fr_tost()` function.")
  }
  
  # Validation for max.print
  if (!is.numeric(max.print) ||
      length(max.print) != 1L ||
      is.na(max.print) ||
      max.print < 1) {
    stop("`max.print` must be a single positive number.")
  }
  
  
  # Basic information needed for summary
  data <- object
  
  redundant <- attr(data, "redundant")
  
  redundancy_str <- if (length(redundant) == 2) {
    paste(redundant, collapse = " and ")
  } else {
    paste0(
      paste(redundant[-length(redundant)], collapse = ", "),
      ", and ",
      redundant[length(redundant)]
    )
  }
  
  alpha <- unique(data$alpha)
  delta_lo <- unique(data$delta_lo)
  delta_hi <- unique(data$delta_hi)
  
  # Delta description, i.e., defined using sd or user specified value
  if (identical(attr(data, "delta_str"), "User specified range")) {
    
    delta_text <- "user specified delta"
    
  } else {
    
    delta_text <- paste0(
      "quantified using ",
      gsub("[^0-9.]", "", attr(data, "delta_str")),
      " times residual SD of model"
    )
  }
  
  # Split data into monoculture, within, and between sets
  sets <- c("Monos", "Within", "Between")
  
  criterion_names <- c(
    Monos = "Monoculture performances",
    Within = "Within-group interactions",
    Between = "Between-group interactions"
  )
  
  criterion_data <- sapply(
    sets,
    function(s) data[data$Set == s, , drop = FALSE],
    simplify = FALSE
  )
  
  # Output width
  con_width <- getOption("width", 80)
  
  cat(strrep("-", con_width), "\n")
  cat(sprintf("\nSummary of functional redundancy test for %s.\n\n",
              redundancy_str))
  cat(strrep("-", con_width), "\n")
  
  # Hypotheses
  cat("\nHypotheses:\n\n")
  cat("  H0: The species are not functionally redundant.\n")
  cat("  HA: The species are functionally redundant.\n")
  
  # Test specification
  cat("\nTest specification:\n\n")
  
  cat(sprintf(
    "  Alpha: %s\n",
    alpha
  ))
  cat(sprintf(
    "  Confidence level: %s%%\n",
    format(100 * (1 - 2 * alpha), trim = TRUE)
  ))
  cat(sprintf(
    "  Equivalence margin: (%s, %s) [%s]\n",
    format(round(delta_lo, 2), trim = TRUE),
    format(round(delta_hi, 2), trim = TRUE),
    delta_text
  ))
  
  # Overall result
  overall_p <- max(data$p_tost, na.rm = TRUE)
  overall_equivalent <- isTRUE(all(data$Equivalent))
  
  cat("\nOverall result:\n\n")
  if (overall_equivalent) {
    cat(sprintf(
     "  Functional redundancy established (p %s).\n",
      format_p(overall_p)
    ))
  } else {
    cat(sprintf(
      "  Functional redundancy not established (p %s).\n",
      format_p(overall_p)
    ))
  }
  
  # Individual criterion-level results
  cat("\nIndividual Criteria:\n")
  
  for (s in sets) {
    d <- criterion_data[[s]]
    # Skip absent criteria, e.g. no Between contrasts
    if (nrow(d) == 0L) {
      next
    }
    criterion_p <- max(d$p_tost, na.rm = TRUE)
    criterion_equivalent <- isTRUE(all(d$Equivalent))
    
    cat(sprintf(
      "\n  %s (p %s): [%s/%s] %s equivalent to 0.\n",
      criterion_names[[s]],
      format_p(criterion_p),
      sum(d$p_tost <= alpha),
      length(d$p_tost),
      if(sum(d$p_tost <= alpha) == 1) "contrast" else "contrasts"
    ))
    
    # Individual contrasts
    n <- nrow(d)
    
    if (n > max.print) {
      cat(sprintf(
        "    Showing %d of %d contrasts. Use `max.print` to show more.\n",
        max.print,
        n
      ))
      d_print <- d[seq_len(max.print), , drop = FALSE]
    } else {
      d_print <- d
    }
    
    # Construct comparison labels
    comparison <- paste(
      d_print$Comm1,
      d_print$Comm2,
      sep = " - "
    )
    
    # Construct display table
    contrast_table <- data.frame(
      Contrast = comparison,
      Estimate = sprintf("%.3f", d_print$Estimate),
      "p-value" = gsub("[^0-9.]", "", sapply(d_print$p_tost, format_p))
    )
    
    # Print table manually so this remains base R
    cat("\n")
    
    # Determine column widths
    widths <- pmax(
      nchar(names(contrast_table)),
      sapply(contrast_table, function(z) max(nchar(z), na.rm = TRUE))
    )
    
    # Header
    header <- paste(
      mapply(
        function(name, width) {
          sprintf(
            paste0("%-", width, "s"),
            name
          )
        },
        names(contrast_table),
        widths
      ),
      collapse = "  "
    )
    
    cat("    ", header, "\n", sep = "")
    
    # Separator
    cat(
      "    ",
      paste(
        mapply(
          function(width) strrep("-", width),
          widths
        ),
        collapse = "  "
      ),
      "\n",
      sep = ""
    )
    
    # Rows
    for (j in seq_len(nrow(contrast_table))) {
      
      row <- paste(
        mapply(
          function(value, width) {
            sprintf(
              paste0("%-", width, "s"),
              value
            )
          },
          as.character(contrast_table[j, ]),
          widths
        ),
        collapse = "  "
      )
      
      cat("    ", row, "\n", sep = "")
    }
    
    if (n > max.print) {
      cat(sprintf(
        "    ... %d additional contrasts not shown.\n",
        n - max.print
      ))
    }
  }
  
  # Footnote
  cat("\n")
  cat(
    "Note:  All tests follow the intersection-union principle and the respective", 
    "overall p-value is the maximum p-value across all lower level constituent tests.\n",
    sep = ""
  )
  
  return(invisible(data))
}


