DI_check_and_fit <- function(fmla, y, block, density, treat, family, data, FG) {
  if(!missing(FG)) FG_ <- FG
  fmla <- paste(fmla)
  fmla <- as.formula(paste0(fmla[2], " ~ ", fmla[3]))
  X_matrix <- model.matrix(fmla, data = data)
  X_check <- DI_matrix_check(X_matrix)
  if(X_check) {
    mod <- glm(formula = fmla, family = family, data = data)
    mod$DIcheck_formula <- fmla
  } else {
    # the lm fit gives NAs for coefficients that cannot be estimated
    fit_to_check <- lm(formula = fmla, data = data)
    coefs_to_check <- coef(fit_to_check)
    cols_to_drop <- which(is.na(coefs_to_check))
    new_model_matrix <- X_matrix[, - cols_to_drop]
    new_data <- data.frame(data, new_model_matrix, check.names = FALSE)
    if(block != "block_zero") {
      original_block_index <- min(grep(block, colnames(new_data)))
    } else {
      original_block_index <- integer(0)
    }
    if(density != "density_zero") {
      original_density_index <- min(grep(density, colnames(new_data)))
    } else {
      original_density_index <- integer(0)
    }
    if(!missing(treat)) {
      original_treat_index <- min(grep(treat, colnames(new_data)))
    } else {
      original_treat_index <- integer(0)
    }
    if(block != "block_zero" | density != "density_zero" | !missing(treat)) {
      new_data <- new_data[, - c(original_block_index,
                                 original_density_index,
                                 original_treat_index)]
    }
    colnames(new_model_matrix) <- gsub("`", "", colnames(new_model_matrix))
    names_new_model_matrix <- paste0("`",colnames(new_model_matrix),"`")
    new_fmla <- as.formula(paste(y, "~", "0+", paste(names_new_model_matrix, collapse = "+")))
    mod <- glm(formula = new_fmla, family = family, data = new_data)
    # RV change
    mod$DIcheck_formula <- fmla
  }
  return(mod)
}

proflik_theta <- function(theta, obj, family, int_terms, prop, DImodel, nSpecies, FGnames) {
  mm <- model.matrix(obj)
  
  if(DImodel %in% c("E","AV")) {
    data_theta_E_AV <- obj$data
    data_theta <- data.frame(mm)
    #data_theta[1:nSpecies] <- data_theta[1:nSpecies]^theta
    new_E_AV <- DI_data_E_AV(prop = prop, data = data_theta_E_AV, theta = theta)
    data_theta_E_AV$E <- new_E_AV$E
    data_theta_E_AV$AV <- new_E_AV$AV
    fitted_model_theta <- glm(formula(obj), family = family, data = data_theta_E_AV)
    #mu_hat <- fitted(fitted_model_theta)
    #sigma_hat <- sqrt(sum(fitted_model_theta$residuals^2)/(fitted_model_theta$df.residual-1))
    #llik <- sum(dnorm(fitted_model_theta$y, mu_hat, sigma_hat, log = TRUE))
    n <- nrow(fitted_model_theta$data)
    p <- length(fitted_model_theta$coef) + 1
    mu_hat <- fitted(fitted_model_theta)
    sigma_hat <- sqrt(sum(fitted_model_theta$residuals^2)/(fitted_model_theta$df.residual) * (n - p)/n)
    llik <- sum(dnorm(fitted_model_theta$y, mu_hat, sigma_hat, log = TRUE))
  } else if(DImodel == "FG") {
    #data_theta_FG <- obj$data
    data_theta_FG <- data.frame(mm, check.names = FALSE)
    colnames(data_theta_FG) <- gsub("`", "", colnames(data_theta_FG))
    # Add original props to calculate interactions
    data_theta_FG <- cbind(data_theta_FG, obj$data[, prop])
    
    #data_theta[1:nSpecies] <- data_theta[1:nSpecies]^theta
    new_FG <- DI_data_FG(prop = prop, FG = FGnames, data = data_theta_FG[, prop], theta = theta)
    FG_ <- new_FG$FG
    ## if we have column names starting with FG_ already
    ## in data_theta_FG, then substitute columns else it's all good
    FG_cols_in_the_data <- grep("FG_", colnames(data_theta_FG))
    if(length(FG_cols_in_the_data) > 0) {
      if(length(FG_cols_in_the_data) != ncol(FG_)) {
        stop("please rename variables beginning with 'FG_'")
      }
      j <- 1
      for(i in FG_cols_in_the_data) {
        data_theta_FG[,i] <- FG_[,j]
        j <- j + 1
      }
    }
    old_formula <- formula(obj)
    new_formula <- paste0(old_formula[2], " ~ ", old_formula[3])
    
    data_theta_FG$y <- obj$y
    names(data_theta_FG)[length(names(data_theta_FG))] <- paste(old_formula[2])
    
    fitted_model_theta <- glm(as.formula(new_formula), family = family, data = data_theta_FG)
    #mu_hat <- fitted(fitted_model_theta)
    #sigma_hat <- sqrt(sum(fitted_model_theta$residuals^2)/(fitted_model_theta$df.residual-1))
    #llik <- sum(dnorm(fitted_model_theta$y, mu_hat, sigma_hat, log = TRUE))
    n <- nrow(fitted_model_theta$data)
    p <- length(fitted_model_theta$coef) + 1
    mu_hat <- fitted(fitted_model_theta)
    sigma_hat <- sqrt(sum(fitted_model_theta$residuals^2)/(fitted_model_theta$df.residual) * (n - p)/n)
    llik <- sum(dnorm(fitted_model_theta$y, mu_hat, sigma_hat, log = TRUE))
  } else if(DImodel == "ADD") {
    #data_theta_ADD <- obj$data
    # Drop existing _add columns to avoid conflicts
    data_theta <- data.frame(mm, check.names = FALSE)
    data_theta <- cbind(data_theta, obj$data[, prop])
    new_ADD <- DI_data_ADD_theta(prop = prop, data = data_theta, theta = theta)
    ADD_cols_in_the_data <- int_terms #grep("_add", colnames(data_theta))
    if(length(ADD_cols_in_the_data) > 0) {
      j <- 1
      for(i in ADD_cols_in_the_data) {
        data_theta[,i] <- new_ADD$ADD_theta[,j]
        j <- j + 1
      }
    }
    old_formula <- formula(obj)
    new_formula <- paste0(old_formula[2], " ~ ", old_formula[3])
    data_theta_ADD <- data_theta
    data_theta_ADD$y <- obj$y
    names(data_theta_ADD)[length(names(data_theta_ADD))] <- paste(old_formula[2])
    names(data_theta_ADD) <- gsub('`','', names(data_theta_ADD))
    fitted_model_theta <- glm(as.formula(new_formula), family = family, data = data_theta_ADD)
    #mu_hat <- fitted(fitted_model_theta)
    #sigma_hat <- sqrt(sum(fitted_model_theta$residuals^2)/(fitted_model_theta$df.residual-1))
    #llik <- sum(dnorm(fitted_model_theta$y, mu_hat, sigma_hat, log = TRUE))
    n <- nrow(fitted_model_theta$data)
    p <- length(fitted_model_theta$coef) + 1
    mu_hat <- fitted(fitted_model_theta)
    sigma_hat <- sqrt(sum(fitted_model_theta$residuals^2)/(fitted_model_theta$df.residual) * (n - p)/n)
    llik <- sum(dnorm(fitted_model_theta$y, mu_hat, sigma_hat, log = TRUE))
  } else {
    #mm[,int_terms] <- nSpecies/(nSpecies - 1) * (nSpecies^2 * mm[,int_terms])^theta
    mm[,int_terms] <- (mm[,int_terms])^theta
    ndata <- obj$data
    ndata$mm <- mm
    fitted_model_theta <- glm(update.formula(formula(obj), . ~ 0 + mm), family = family, data = ndata)
    #mu_hat <- fitted(fitted_model_theta)
    #sigma_hat <- sqrt(sum(fitted_model_theta$residuals^2)/(fitted_model_theta$df.residual-1))
    #llik <- sum(dnorm(fitted_model_theta$y, mu_hat, sigma_hat, log = TRUE))
    n <- nrow(fitted_model_theta$data)
    p <- length(fitted_model_theta$coef) + 1
    mu_hat <- fitted(fitted_model_theta)
    sigma_hat <- sqrt(sum(fitted_model_theta$residuals^2)/(fitted_model_theta$df.residual) * (n - p)/n)
    llik <- sum(dnorm(fitted_model_theta$y, mu_hat, sigma_hat, log = TRUE))
  }
  return(llik)
}

get_int_terms_FULL <- function(mm_names, prop_names) {
  # removing ` characters
  mm_names <- gsub("`", "", mm_names)
  # finding all terms that include one and only one ":" operator
  all_pairwise_ints <- which(lengths(regmatches(mm_names,
                                                gregexpr(":", mm_names))) == 1)
  # matching species names to each term in the model matrix
  prop_match <- lapply(strsplit(mm_names, ":"),
                       function(x) {
                         x %in% prop_names
                       })
  # finding all terms that include a pairwise interaction between species
  prop_pairwise <- which(unlist(lapply(prop_match, sum)) == 2)
  # the answer is the intersection between all_pairwise_ints and prop_pairwise,
  # i.e., a pairwise interaction between species and nothing else
  int_terms <- intersect(all_pairwise_ints, prop_pairwise)
  return(int_terms)
}

DI_theta <- function(obj, DImodel, FGnames, prop, nSpecies, family) {
  if(missing(FGnames)) {
    FGnames <- NULL
  }
  mm <- model.matrix(obj)
  int_terms <- switch(EXPR = DImodel,
                      "AV" = grep("AV", colnames(mm)),
                      "E" = grep("E", colnames(mm)),
                      "FG" = grep("FG", colnames(mm)),
                      "ADD" = grep("_add", colnames(mm)),
                      "FULL" = get_int_terms_FULL(mm_names = colnames(mm),
                                                  prop_names = prop))
  #options(warn = -1)
  upper_boundary <- 1.5
  theta_info <- get_theta_info(upper_boundary = upper_boundary, DImodel = DImodel,
                               obj = obj, prop = prop, 
                               family = family, int_terms = int_terms, 
                               nSpecies = nSpecies, FGnames = FGnames)
  #options(warn = 0)
  theta_hat <- theta_info$theta_hat
  profile_loglik <- theta_info$profile_loglik
  if((upper_boundary - theta_hat) < .01) {
    warning("Theta has reached the upper boundary, this may indicate lack of convergence and/or a problem with non-significant diversity effect.")
  }
  
  if(DImodel %in% c("E","AV")) {
    data_theta_E_AV <- obj$data
    data_theta <- data.frame(mm)
    #data_theta[1:nSpecies] <- data_theta[1:nSpecies]^theta_hat
    new_E_AV <- DI_data_E_AV(prop = prop, data = data_theta_E_AV, theta = theta_hat)
    data_theta_E_AV$E <- new_E_AV$E
    data_theta_E_AV$AV <- new_E_AV$AV
    mod_theta <- glm(formula(obj), family = family, data = data_theta_E_AV)
  } else if(DImodel == "FG") {
    #data_theta_FG <- obj$data
    data_theta_FG <- data.frame(mm, check.names = FALSE)
    colnames(data_theta_FG) <- gsub("`", "", colnames(data_theta_FG))
    data_theta_FG <- cbind(data_theta_FG, obj$data[, prop])
    #data_theta[1:nSpecies] <- data_theta[1:nSpecies]^theta_hat
    new_FG <- DI_data_FG(prop = prop, FG = FGnames, theta = theta_hat, data = data_theta_FG)
    FG_ <- new_FG$FG
    ## if we have column names starting with FG_ already
    ## in data_theta_FG, then substitute columns, else it's all good
    FG_cols_in_the_data <- grep("FG_", colnames(data_theta_FG))
    if(length(FG_cols_in_the_data) > 0) {
      j <- 1
      for(i in FG_cols_in_the_data) {
        data_theta_FG[,i] <- FG_[,j]
        j <- j + 1
      }
    }
    old_formula <- formula(obj)
    new_formula <- paste0(old_formula[2], " ~ ", old_formula[3])
    data_theta_FG$y <- obj$y
    names(data_theta_FG)[length(names(data_theta_FG))] <- paste(old_formula[2])
    mod_theta <- glm(as.formula(new_formula), family = family, data = data_theta_FG)
  } else if(DImodel == "ADD") {
    #data_theta_ADD <- obj$data
    data_theta <- data.frame(mm, check.names = FALSE)
    data_theta <- cbind(data_theta, obj$data[, prop])
    new_ADD <- DI_data_ADD_theta(prop = prop, data = data_theta, theta = theta_hat)
    ADD_cols_in_the_data <- int_terms #grep("_add", colnames(data_theta))
    if(length(ADD_cols_in_the_data) > 0) {
      j <- 1
      for(i in ADD_cols_in_the_data) {
        data_theta[,i] <- new_ADD$ADD_theta[,j]
        j <- j + 1
      }
    }
    old_formula <- formula(obj)
    new_formula <- paste0(old_formula[2], " ~ ", old_formula[3])
    data_theta_ADD <- data_theta
    data_theta_ADD$y <- obj$y
    names(data_theta_ADD)[length(names(data_theta_ADD))] <- paste(old_formula[2])
    names(data_theta_ADD) <- gsub('`','', names(data_theta_ADD))
    mod_theta <- glm(as.formula(new_formula), family = family, data = data_theta_ADD)
  } else {
    #mm[,int_terms] <- nSpecies/(nSpecies - 1) * (nSpecies^2 * mm[,int_terms])^theta_hat
    mm[,int_terms] <- (mm[,int_terms])^theta_hat
    colnames(mm) <- gsub("`", "", colnames(mm))
    names_mm <- paste0("`", colnames(mm)[int_terms], "`")
    resp_name <- paste(formula(obj))[2]
    ndata <- data.frame(obj$data[,resp_name], mm, check.names = FALSE)
    names(ndata)[1] <- resp_name
    colnames(mm) <- paste0("`",colnames(mm),"`")
    new_formula_theta <- as.formula(paste(resp_name, "~", "0+",
                                          paste(colnames(mm)[-int_terms], collapse = "+"), "+",
                                          paste(names_mm, collapse = "+")))
    mod_theta <- glm(formula = new_formula_theta,
                     family = family, data = ndata)
  }
  
  mod_theta$coefficients <- c(mod_theta$coefficients, "theta" = theta_hat)
  mod_theta$df.residual <- mod_theta$df.residual - 1
  mod_theta$profile_loglik <- profile_loglik
  mod_theta$aic <- AIC2(mod_theta)
  return(mod_theta)
}

get_theta_info <- function(upper_boundary, DImodel, obj, prop, family, int_terms,
                           nSpecies, FGnames) {
  optimum <- optimize(proflik_theta, interval = c(0.00001, upper_boundary), 
                      maximum = TRUE, DImodel = DImodel, obj = obj, prop = prop, family = family, 
                      int_terms = int_terms, nSpecies = nSpecies, FGnames = FGnames)
  theta_hat <- optimum$maximum
  theta_grid <- seq(0.00001, upper_boundary + 1, length = 100)
  proflik_theta_vec <- Vectorize(proflik_theta, "theta")
  profile_loglik <- proflik_theta_vec(theta = theta_grid, obj = obj, prop = prop,
                                      family = family, int_terms = int_terms, DImodel = DImodel, 
                                      nSpecies = nSpecies, FGnames = FGnames)
  
  # Adding theta_hat in the profile likelihood grid used for searching in the CI.
  # This will ensure that the CI is always calculated on the estimate of theta
  profile_loglik = data.frame(prof = c(profile_loglik, optimum$objective), 
                              grid = c(theta_grid,theta_hat))
  profile_loglik <- profile_loglik[order(profile_loglik$grid),]
  return(list("theta_hat" = theta_hat,
              "profile_loglik" = profile_loglik))
}

theta_CI <- function(obj, conf = .95, n = 100) {
  threshold <- max(obj$profile_loglik$prof) - qchisq(conf, 1)/2
  if(threshold < min(obj$profile_loglik$prof) | threshold > max(obj$profile_loglik$prof)) {
    stop("CI cannot be computed. This is because the profile log-likelihood function is flat or displays unusual behaviour at the interval theta = (0.01, 2.5).") 
  }
  CI_finder <- approxfun(x = obj$profile_loglik$grid,
                         y = obj$profile_loglik$prof - threshold)
  
  CI <- rootSolve::uniroot.all(CI_finder, interval = range(obj$profile_loglik$grid), n = n)
  
  # Increasing n if convergence fails
  if(length(CI) == 0){
    CI <- rootSolve::uniroot.all(CI_finder, interval = range(obj$profile_loglik$grid), n = 100000)
  }
  
  alpha <- 1 - conf
  
  # Throw error if only one root can be returned
  if(length(CI)!=2){
    stop(paste0('There are convergence problems and a ', conf*100, '% CI couldn\'t be computed. Try using a lower confidence level.'))
  }
  
  names(CI) <- c("lower","upper")
  return(CI)
}

namesub_DI <- Vectorize(function(name) {
  thename <- switch(name,
                    "CUSTOM" = "Custom DI model",
                    "STR" = "Structural 'STR' DImodel",
                    "ID" = "Species identity 'ID' DImodel",
                    "AV" = "Average interactions 'AV' DImodel",
                    "E" = "Evenness 'E' DImodel",
                    "ADD" = 
                      "Additive species contributions to interactions 'ADD' DImodel",
                    "FG" = "Functional group effects 'FG' DImodel",
                    "FULL" = "Separate pairwise interactions 'FULL' DImodel",
                    "STR_treat" = "Structural 'STR' DImodel with treatment covariate",
                    "ID_treat" = "Species identity 'ID' DImodel with treatment covariate",
                    "AV_treat" = "Average interactions 'AV' DImodel with treatment covariate",
                    "E_treat" = "Evenness 'E' DImodel with treatment covariate",
                    "ADD_treat" = 
                      "Additive species contributions to interactions 'ADD' DImodel with treatment covariate",
                    "FG_treat" = "Functional group effects 'FG' DImodel with treatment covariate",
                    "FULL_treat" = 
                      "Separate pairwise interactions 'FULL' DImodel with treatment covariate",
                    stop("not yet implemented"))
  return(thename)
}, "name")

get_community <- function(prop, data) {
  Pind <- get_P_indices(prop = prop, data = data)$Pind
  comms <- data[, Pind]
  n_obs <- nrow(comms)
  unique_comms <- unique(comms)
  n_comms <- nrow(unique_comms)
  ## function that identifies whether all values are the same in a community
  identifier_fun <- function(x) {
    apply(comms, 1, function(y) all(y == x))
  }
  ## matrix indicating which community is in which row
  id_matrix <- apply(unique_comms, 1, identifier_fun)
  ## community id vector
  comm_id <- apply(id_matrix, 1, which)
  ## transforming into a factor
  community_factor <- as.factor(comm_id)
  ## message and return
  message("'community' is a factor with ", n_comms, " levels, one for each unique set of proportions.")
  return(community_factor)
}

DI_compare <- function(model, ...) {
  theta_flag <- attr(model, "theta_flag")
  og_data <- model$original_data
  prop <- attr(model, "prop")
  if(theta_flag) { # retrieve theta value
    theta <- coef(model)["theta"]    # retrieve theta value
  } else {                           # retrieve theta value
    theta <- 1                       # retrieve theta value
  }                                  # retrieve theta value
  ref_model <- DI_reference(model = model, prop = prop, 
                            data = og_data, theta = theta) # add theta argument
  print(anova(model, ref_model, ...))
  return(invisible(ref_model)) # return ref_model instead of function DI_reference
}

DI_reference <- function(model, prop, data, theta = 1) { # add theta argument
  community <- get_community(prop = prop, data = data)
  new_data <- data
  new_data$community <- community
  # delete line replacing call with DIcall
  ref_model <- update_DI(model, # use update_DI function, no need to replace call
                         extra_formula = ~ community,
                         estimate_theta = FALSE, # turn off theta estimation
                         theta = theta,          # use estimated value for original model
                         data = new_data)
  return(ref_model)
}

anova.DI <- function(object, ...) {
  input <- as.list(match.call())
  #print(input)
  input <- input[-1]
  #print(input)
  if(length(which(names(input) == "test")) > 0) input <- input[- which(names(input) == "test")]
  if(length(input) == 1) {
    stop("anova method not yet implemented for single DI model objects. You can only use the anova function to compare multiple nested DI models.")
  } else {
    anovaDIglm(object, ...) 
  }
}

anovaDIglm <- function (object, ..., dispersion = NULL, test = NULL) {
  dotargs <- list(...)
  named <- if (is.null(names(dotargs))) 
    rep_len(FALSE, length(dotargs))
  else (names(dotargs) != "")
  if (any(named)) 
    warning("the following arguments to 'anova.glm' are invalid and dropped: ", 
            paste(deparse(dotargs[named]), collapse = ", "))
  dotargs <- dotargs[!named]
  is.glm <- vapply(dotargs, function(x) inherits(x, "glm"), 
                   NA)
  dotargs <- dotargs[is.glm]
  if (length(dotargs)) 
    return(anova_glmlist(c(list(object), dotargs), dispersion = dispersion, 
                         test = test))
  doscore <- !is.null(test) && test == "Rao"
  varlist <- attr(object$terms, "variables")
  x <- if (n <- match("x", names(object), 0L)) 
    object[[n]]
  else model.matrix(object)
  varseq <- attr(x, "assign")
  nvars <- max(0, varseq)
  resdev <- resdf <- NULL
  if (doscore) {
    score <- numeric(nvars)
    method <- object$method
    y <- object$y
    fit <- eval(call(if (is.function(method)) "method" else method, 
                     x = x[, varseq == 0, drop = FALSE], y = y, weights = object$prior.weights, 
                     start = object$start, offset = object$offset, family = object$family, 
                     control = object$control))
    r <- fit$residuals
    w <- fit$weights
    icpt <- attr(object$terms, "intercept")
  }
  if (nvars > 1 || doscore) {
    method <- object$method
    y <- object$y
    if (is.null(y)) {
      mu.eta <- object$family$mu.eta
      eta <- object$linear.predictors
      y <- object$fitted.values + object$residuals * mu.eta(eta)
    }
    for (i in seq_len(max(nvars - 1L, 0))) {
      fit <- eval(call(if (is.function(method)) "method" else method, 
                       x = x[, varseq <= i, drop = FALSE], y = y, weights = object$prior.weights, 
                       start = object$start, offset = object$offset, 
                       family = object$family, control = object$control))
      if (doscore) {
        zz <- eval(call(if (is.function(method)) "method" else method, 
                        x = x[, varseq <= i, drop = FALSE], y = r, 
                        weights = w, intercept = icpt))
        score[i] <- zz$null.deviance - zz$deviance
        r <- fit$residuals
        w <- fit$weights
      }
      resdev <- c(resdev, fit$deviance)
      resdf <- c(resdf, fit$df.residual)
    }
    if (doscore) {
      zz <- eval(call(if (is.function(method)) "method" else method, 
                      x = x, y = r, weights = w, intercept = icpt))
      score[nvars] <- zz$null.deviance - zz$deviance
    }
  }
  resdf <- c(object$df.null, resdf, object$df.residual)
  resdev <- c(object$null.deviance, resdev, object$deviance)
  table <- data.frame(c(NA, -diff(resdf)), c(NA, pmax(0, -diff(resdev))), 
                      resdf, resdev)
  tl <- attr(object$terms, "term.labels")
  if (length(tl) == 0L) 
    table <- table[1, , drop = FALSE]
  dimnames(table) <- list(c("NULL", tl), c("Df", "Deviance", 
                                           "Resid. Df", "Resid. Dev"))
  if (doscore) 
    table <- cbind(table, Rao = c(NA, score))
  title <- paste0("Analysis of Deviance Table", "\n\nModel: ", 
                  object$family$family, ", link: ", object$family$link, 
                  "\n\nResponse: ", as.character(varlist[-1L])[1L], "\n\nTerms added sequentially (first to last)\n\n")
  df.dispersion <- Inf
  if (is.null(dispersion)) {
    dispersion <- summary(object, dispersion = dispersion)$dispersion
    df.dispersion <- if (dispersion == 1) 
      Inf
    else object$df.residual
  }
  if (!is.null(test)) {
    if (test == "F" && df.dispersion == Inf) {
      fam <- object$family$family
      if (fam == "binomial" || fam == "poisson") 
        warning(gettextf("using an F test with a '%s' family is inappropriate", 
                         fam), domain = NA)
      else warning("using an F test with a fixed dispersion is inappropriate")
    }
    table <- stat.anova(table = table, test = test, scale = dispersion, 
                        df.scale = df.dispersion, n = NROW(x))
  }
  structure(table, heading = title, class = c("anova", "data.frame"))
}

anova_glmlist <- function (object, ..., dispersion = NULL, test = NULL) {
  doscore <- !is.null(test) && test == "Rao"
  responses <- as.character(lapply(object, function(x) {
    deparse(formula(x)[[2L]])
  }))
  sameresp <- responses == responses[1L]
  if (!all(sameresp)) {
    object <- object[sameresp]
    warning(gettextf("models with response %s removed because response differs from model 1", 
                     sQuote(deparse(responses[!sameresp]))), domain = NA)
  }
  ns <- sapply(object, function(x) length(x$residuals))
  if (any(ns != ns[1L])) 
    stop("models were not all fitted to the same size of dataset")
  nmodels <- length(object)
  if (nmodels == 1) 
    return(anovaDIglm(object[[1L]], dispersion = dispersion, 
                     test = test))
  resdf <- as.numeric(lapply(object, function(x) x$df.residual))
  resdev <- as.numeric(lapply(object, function(x) x$deviance))

  if (doscore) {
    score <- numeric(nmodels)
    score[1] <- NA
    df <- -diff(resdf)
    for (i in seq_len(nmodels - 1)) {
      m1 <- if (df[i] > 0) 
        object[[i]]
      else object[[i + 1]]
      m2 <- if (df[i] > 0) 
        object[[i + 1]]
      else object[[i]]
      r <- m1$residuals
      w <- m1$weights
      method <- m2$method
      icpt <- attr(m1$terms, "intercept")
      zz <- eval(call(if (is.function(method)) "method" else method, 
                      x = model.matrix(m2), y = r, weights = w, intercept = icpt))
      score[i + 1] <- zz$null.deviance - zz$deviance
      if (df[i] < 0) 
        score[i + 1] <- -score[i + 1]
    }
  }
  table <- data.frame(resdf, resdev, c(NA, -diff(resdf)), 
                      c(NA, -diff(resdev)))
  variables <- lapply(object, function(x) paste(deparse(formula(x)), 
                                                collapse = "\n"))
  dimnames(table) <- list(1L:nmodels, c("Resid. Df", "Resid. Dev", 
                                        "Df", "Deviance"))
  if (doscore) 
    table <- cbind(table, Rao = score)
  title <- "Analysis of Deviance Table\n"
  topnote <- paste0("Model ", format(1L:nmodels), ": ", variables, 
                    collapse = "\n")
  if (!is.null(test)) {
    bigmodel <- object[[order(resdf)[1L]]]
    dispersion <- summary(bigmodel, dispersion = dispersion)$dispersion
    df.dispersion <- if (dispersion == 1) 
      Inf
    else min(resdf)
    if (test == "F" && df.dispersion == Inf) {
      fam <- bigmodel$family$family
      if (fam == "binomial" || fam == "poisson") 
        warning(gettextf("using an F test with a '%s' family is inappropriate", 
                         fam), domain = NA, call. = FALSE)
      else warning("using an F test with a fixed dispersion is inappropriate")
    }
    table <- stat.anova(table = table, test = test, scale = dispersion, 
                        df.scale = df.dispersion, n = length(bigmodel$residuals))
  }
  structure(table, heading = c(title, topnote), class = c("anova", 
                                                          "data.frame"))
}

update_DI <- function(object, ...) {
  object$call <- object$DIcall
  update.default(object, ...)
}