DI_data_prepare <- function(y, block, density, prop, treat, ID = NULL, FG = NULL, data, theta = 1) {
  if(any(is.na(data))) {
    stop("The dataset contains missing values. Please remove them prior to analysis.") 
  }
  if(missing(y)) {
    y <- rep(0, nrow(data))
    warning("y was not supplied, so a vector of zeros was created for the response\n")
  }
  if(missing(block)) block <- NA
  if(missing(density)) density <- NA
  if(missing(treat)) treat <- NA
  # variables directly obtained from dataset
  if(is.na(block)) {
    block <- "block_zero"
    data$block_zero <- 0
  } else {
    if(!is.character(block)) block <- names(data)[block] # variable block (name)
    ## warning for continuous block variables
    if(!is.factor(data[,block])) {
      #data[,block] <- as.factor(data[,block])
      warning("'", block, "' has been declared as a block variable, but is not a factor.")#,
              #" It has been converted to a factor.")
    }
  }
  if(is.na(density)) {
    density <- "density_zero"
    data$density_zero <- 0
  } else {
    if(!is.character(density)) density <- names(data)[density] # variable density (name)
    ## warning for continuous density variables
    if(!is.factor(data[,density])) {
      #data[,density] <- as.factor(data[,density])
      warning("'", density, "' has been declared as a density variable, but is not a factor.")#,
      #" It has been converted to a factor.")
    }
  }
  if(is.na(treat)) {
    treat <- "treat_zero"
    data$treat_zero <- 0
  }
  
  if(!is.character(y)) y <- names(data)[y] # response variable (name)
  
  ## getting the indices for the columns corresponding to the species proportions (Pind)
  Pind_and_prop <- get_P_indices(prop = prop, data = data)
  Pind <- Pind_and_prop$Pind
  
  ## initial proportion checks - looking for rows where all species proportions = 0
  prop_zero_rows <- get_zero_prop(Pind = Pind, data = data)
  if(length(prop_zero_rows) > 0) {
    warning("One or more rows in your dataset have ALL proportions equal to zero. They are included in the analysis, but if these rows are in your dataset in error this will affect the results.")
    data_check <- data[- prop_zero_rows,]
  } else {
    data_check <- data
  }
  
  ## checking if the Pi's sum to 1
  prop_check <- DI_prop_check(Pind = Pind, data = data_check)
  if(prop_check == "OOB") {
    stop("One or more rows have species proportions with values less than 0 or greater than 1. This must be corrected prior to analysis.\n")
  }
  if(prop_check == "error") {
    stop("One or more rows have species proportions that do not sum to 1. This must be corrected prior to analysis.\n")
  }
  if(prop_check == "minor") {
    Pi_sums <- apply(data[,Pind], 1, sum)
    Pi_sums <- ifelse(Pi_sums == 0, 1, Pi_sums)
    data[,Pind] <- data[,Pind]/Pi_sums
    warning("One or more rows have species proportions that sum to approximately 1, but not exactly 1. This is typically a rounding issue, and has been corrected internally prior to analysis.\n")
  }
  
  ## species proportions column names (prop)
  prop <- Pind_and_prop$prop

  if(!is.character(treat)) treat <- names(data)[treat] # environmental covariate (name)
  ## calculating the E and AV variables
  E_AV <- DI_data_E_AV(prop = prop, data = data, theta = theta)
  newdata <- data
  newdata$AV <- E_AV$AV
  newdata$E <- E_AV$E
  even_flag <- E_AV$even_flag
  ## calculating the P*(1-P) (P_add) variables
  ADD <- DI_data_ADD(prop = prop, data = data, theta = theta)
  P_int_flag <- ADD$P_int_flag
  if(P_int_flag) {
    ADD$ADD <- 0
  }
  newdata <- data.frame(newdata, ADD$ADD)
  ## calculating FG variables
  if(!is.null(FG)) {
    FG <- DI_data_FG(prop = prop, FG = FG, data = data, theta = theta)$FG
    FG_names <- colnames(FG)
  } else {
    FG_names <- NULL
  }
  
  # Calculating the FULL variables
  FULL <- DI_data_fullpairwise(prop = prop, data = data[, prop], theta = theta)
  # Add variabes to data
  newdata <- cbind(newdata, FULL)
  
  ###### RV change########
  # Grouping ID effects
  if(is.null(ID)){
    ID <- paste0(prop, "_ID")
  } 
  ID_name_check(ID = ID, prop = prop, FG = FG_names)
  grouped_ID <- group_IDs(data = data, prop = prop, ID = ID)
  newdata <- cbind(newdata, grouped_ID)
  #######################
  
  ## return object
  return(list("newdata" = newdata, "y" = y, "block" = block, density = density,
              "prop" = prop, "ID" = unique(ID), ID_map = ID, "treat" = treat, "FG" = FG,
              "P_int_flag" = P_int_flag, "even_flag" = even_flag, "nSpecies" = length(prop)))
}

get_P_indices <- function(prop, data) {
  if(!is.character(prop)) {
    Pind <- prop # indices for the columns with species proportions P_i
    prop <- names(data[, prop]) # species proportions P_i (names)
  } else {
    vec_grep <- Vectorize(grep, "pattern")
    Pind <- as.numeric(vec_grep(paste0("^", prop, "$"),#paste("\\<", prop, "\\>", sep = ""),
                                names(data)))
  }
  return(list(Pind = Pind, prop = prop))
}

DI_data_E_AV_internal <- function(prop, data) {
  if(!is.character(prop)) {
    prop <- get_P_indices(prop = prop, data = data)$prop
  }
  nSpecies <- length(prop) # number of species
  if(nSpecies <= 1) stop("must have at least 2 species to fit DI models")
  nComb <- choose(nSpecies, 2) # number of pairwise combinations
  pairwise_fmla <- as.formula(paste("~", "0+", "(",
                                    paste(prop, collapse = "+"), ")^2"))
  normE <- 2 * nSpecies/(nSpecies - 1)
  Ematrix <- model.matrix(pairwise_fmla, data = data)
  # evenness and "AV" variable
  if(nSpecies > 2) {
    AV <- rowSums(Ematrix[,(nSpecies+1):ncol(Ematrix)])
    E <- normE * AV
    even_flag <- FALSE
    return(list("AV" = AV, "E" = E, "even_flag" = even_flag))
  } else {
    even_flag <- TRUE
    return(list("even_flag" = even_flag))
  }
}

DI_data_E_AV <- function(prop, data, theta = 1) {
  ########## RV change ###############
  if (inherits(data, 'tbl')){
    data <- as.data.frame(data)
  }
  
  ####################################
  if(theta != 1) {
    data[,prop] <- data[,prop]^theta 
  }
  result <- DI_data_E_AV_internal(prop = prop, data = data)
  return(result)
}

DI_data_ADD_internal <- function(prop, data) {
  Pind <- get_P_indices(prop = prop, data = data)$Pind
  nSpecies <- length(prop)
  if(nSpecies > 3) {
    ADD_vars <- as.data.frame(apply(data[, Pind], 2, function(x) x*(1-x)))
    prop_names <- names(data)[Pind]
    names(ADD_vars) <- paste(prop_names, "_add", sep = "")
    P_int_flag <- FALSE
  } else {
    ADD_vars <- NULL
    P_int_flag <- TRUE
  }
  return(list("ADD" = ADD_vars, "P_int_flag" = P_int_flag))
}

DI_data_ADD_theta <- function(prop, data, theta) {
  Pind <- get_P_indices(prop = prop, data = data)$Pind
  nSpecies <- length(prop)
  if(nSpecies > 3) {
    P_matrix <- data[,Pind]
    Pi_theta <- P_matrix^theta
    ADD_vars_theta <- matrix(NA, ncol = nSpecies, nrow = nrow(data))
    for(i in 1:nSpecies) {
      sum_Pj_theta <- apply(Pi_theta[,-i], 1, sum) 
      ADD_vars_theta[,i] <- Pi_theta[,i] * sum_Pj_theta
    }
    ADD_vars_theta <- as.data.frame(ADD_vars_theta)
    prop_names <- names(data)[Pind]
    names(ADD_vars_theta) <- paste(prop_names, "_add", sep = "")
    P_int_flag <- FALSE
  } else {
    ADD_vars_theta <- NULL
    P_int_flag <- TRUE
  }
  return(list("ADD_theta" = ADD_vars_theta, "P_int_flag" = P_int_flag))
}

DI_data_ADD <- function(prop, data, theta = 1) {
  ########## RV change ###############
  if (inherits(data, 'tbl')){
    data <- as.data.frame(data)
  }
  add_name_check(data = data)
  ####################################
  if(theta == 1) {
    result <- DI_data_ADD_internal(prop = prop, data = data)
  } else {
    result <- DI_data_ADD_theta(prop = prop, data = data, theta = theta) 
  }
  return(result)
}

DI_data_FG_internal <- function(prop, FG, data) {
  # name check
  FG_name_check(FG = FG)
  # number of functional groups
  nfg <- length(unique(FG))
  # n for checking at the end
  n_check <- choose(nfg, 2) + nfg
  
  Pind <- get_P_indices(prop = prop, data = data)$Pind
  if(any(!is.character(FG))) stop("FG argument takes character strings with functional",
                                  " group names referring to each species, in order")
  fg_index <- FG
  fg_index_names <- levels(as.factor(fg_index))
  testdata <- data[,Pind]
  prop <- colnames(testdata)
  prop <- paste(prop, fg_index, sep = "_infg_")
  testdata2 <- testdata
  colnames(testdata2) <- prop
  nSpecies <- length(prop)
  pairwise_fmla <- as.formula(paste("~", "0+", "(",
                                    paste(prop, collapse = "+"), ")^2"))
  Pmatrix <- model.matrix(pairwise_fmla, data = testdata2)
  FGmatrix_raw <- Pmatrix[,(nSpecies+1):ncol(Pmatrix)]
  FGrawnames <- colnames(FGmatrix_raw)
  ## matching between FG effects
  #nfg <- length(fg_index_names)
  FGmatch <- list()
  for(i in 1:nfg) {
    FGmatch[[i]] <- grep(fg_index_names[i], FGrawnames)
  }
  FGnames <- rep(NA, length(FGrawnames))
  for(i in 1:nfg) {
    for(j in 1:nfg) {
      if(i < j) {
        FGnames[FGmatch[[i]][FGmatch[[i]] %in% FGmatch[[j]]]] <- 
          paste("bfg_", fg_index_names[i], ":", fg_index_names[j], sep = "")
      }
    }
  }
  ## matching within FG effects
  FGmatch2 <- list()
  for(i in 1:nfg) {
    FGmatch2[[i]] <- grep(paste(fg_index_names[i], ":", sep = ""), FGrawnames)
  }
  for(i in 1:nfg) {
    helper_index <- FGmatch[[i]][FGmatch[[i]] %in% FGmatch2[[i]]]
    FGnames[helper_index[is.na(FGnames[helper_index])]] <- 
      paste("wfg_", fg_index_names[i], sep = "")
  }
  ## summing variables belonging to same FG effect
  FGeffects <- levels(as.factor(FGnames))
  FG <- matrix(NA, nrow = nrow(FGmatrix_raw), ncol = length(FGeffects))
  for(i in 1:length(FGeffects)) {
    FG[,i] <- rowSums(as.matrix(FGmatrix_raw[,FGnames == FGeffects[i]]))
  }
  ## names of FG variables
  FGeffects <- gsub(":","_", FGeffects)
  colnames(FG) <- FGeffects

  ## final checks
  if(any(table(fg_index) < 2)) {
    return(list("FG" = FG))
  } else if(n_check != ncol(FG)) {
    stop("Expected ", n_check, " terms, but have ", ncol(FG),
         ". Please give your functional groups a different name.",
         " One or more of the options are being used internally.",
         " Perhaps use upper-case names?")
  }
  
  return(list("FG" = FG))
}

DI_data_FG <- function(prop, FG, data, theta = 1) {
  ########## RV change ###############
  if (inherits(data, 'tbl')){
    data <- as.data.frame(data)
  }
  
  ####################################
  if(theta != 1) {
    data[,prop] <- data[,prop]^theta 
  }
  result <- DI_data_FG_internal(prop = prop, FG = FG, data = data)
  return(result)
}

DI_data_fullpairwise <- function(prop, data, theta = 1) {
  ########## RV change ###############
  if (inherits(data, 'tbl')){
    data <- as.data.frame(data)
  }
  
  ####################################
  prop <- get_P_indices(prop = prop, data = data)$prop
  fmla <- as.formula(paste("~ 0 + ", "(", paste(prop, collapse = "+"), ")^2"))
  all_pairwise <- model.matrix(fmla, data = data)
  n_species <- length(prop)
  obj <- all_pairwise[,-c(1:n_species)]^theta
  
  # Special case for 2 species system
  if(n_species == 2){
    int <- paste(prop, collapse = ':')
    obj <- matrix(obj, ncol = 1, dimnames = list(NULL, c(int)))
    names(obj) <- int
  }
  return(obj)
}

DI_prop_check <- function(Pind, data) {
  props <- data[,Pind]
  # prop values are out of bounds, i.e. > 1 or < 0
  if(any(props > 1 | props < 0)){
    return("OOB") 
  }
  pi_sums <- apply(props, 1, sum)
  if(any(pi_sums > 1.0001) | any(pi_sums < .9999)) {
    return("error")
  } else if(any(pi_sums < 1 & pi_sums > .9999) | any(pi_sums > 1 & pi_sums < 1.0001)) {
    return("minor") 
  } else return("ok")
}

get_zero_prop <- function(Pind, data) {
  props <- data[,Pind]
  pi_zero <- apply(props, 1, function(x) all(x == 0))
  return(which(pi_zero))
}

FG_name_check <- function(FG) {
  cond1 <- length(grep(":", FG)) > 0
  cond2 <- any(FG == "_")
  cond3 <- any(FG == "i")
  cond4 <- any(FG == "n")
  cond5 <- any(FG == "f")
  cond6 <- any(FG == "g")
  cond7 <- any(FG == "_i")
  cond8 <- any(FG == "in")
  cond9 <- any(FG == "nf")
  cond10 <- any(FG == "fg")
  cond11 <- any(FG == "g_")
  cond12 <- any(FG == "_in")
  cond13 <- any(FG == "inf")
  cond14 <- any(FG == "nfg")
  cond15 <- any(FG == "fg_")
  cond16 <- any(FG == "_inf")
  cond17 <- any(FG == "infg")
  cond18 <- any(FG == "nfg_")
  cond19 <- any(FG == "_infg")
  cond20 <- any(FG == "infg_")
  cond21 <- any(FG == "_infg_")
  
  if(cond1 | cond2 | cond3 | cond4 | cond5 | cond6 | cond7 |
     cond8 | cond9 | cond10 | cond11 | cond12 | cond13 | cond14 |
     cond15 | cond16 | cond17 | cond18 | cond19 | cond20 | cond21) {
    stop("Please give your functional groups a different name.",
         " Names should not include colons (':'), or any single or multiple",
         " character combination of the expression '_infg_'.",
         " This expression is reserved for computing functional groups internally.")
  }
}

############# RV change #################
add_name_check <- function(data){
  if (any(grepl('^.*_add', colnames(data)))){
    stop('Certain column names cause internal conflicts when calculating additive interactions. Please rename any columns containing the string \'_add\'.')
  }
}

ID_name_check <- function(ID, prop, FG){
  cond1 <- length(grep(":", ID)) > 0
  cond2 <- any(ID == "_")
  cond3 <- any(ID == "i")
  cond4 <- any(ID == "n")
  cond5 <- any(ID == "f")
  cond6 <- any(ID == "g")
  cond7 <- any(ID == "_i")
  cond8 <- any(ID == "in")
  cond9 <- any(ID == "nf")
  cond10 <- any(ID == "fg")
  cond11 <- any(ID == "g_")
  cond12 <- any(ID == "_in")
  cond13 <- any(ID == "inf")
  cond14 <- any(ID == "nfg")
  cond15 <- any(ID == "fg_")
  cond16 <- any(ID == "_inf")
  cond17 <- any(ID == "infg")
  cond18 <- any(ID == "nfg_")
  cond19 <- any(ID == "_infg")
  cond20 <- any(ID == "infg_")
  cond21 <- any(ID == "_infg_")
  cond22 <- length(grep(pattern = "^add", x = ID, ignore.case = TRUE)) > 0
  cond23 <- length(grep(pattern = "^full", x = ID, ignore.case = TRUE)) > 0
  cond24 <- length(grep(pattern = "^fg", x = ID, ignore.case = TRUE)) > 0
  cond25 <- any(ID == "E")
  cond26 <- any(ID == "AV")
  if(cond1 | cond2 | cond3 | cond4 | cond5 | cond6 | cond7 |
     cond8 | cond9 | cond10 | cond11 | cond12 | cond13 | cond14 |
     cond15 | cond16 | cond17 | cond18 | cond19 | cond20 | cond21 | 
     cond22 | cond23 | cond24| cond25 | cond26) {
    stop("Please give your IDs a different name.",
         " Names should not include colons (':'), or any single or multiple",
         " character combination of the expressions '_infg_', 'ADD', 'FG', 'FULL', 'E', and 'AV'.",
         " These expressions are reserved for internal computations.")
  }
  if(length(ID)!=length(prop)){
    stop("'ID' should be a vector of same length as prop vector specifying the grouping structure for the ID effects")
  }
  if(any(ID %in% FG)){
    stop("'ID' cannot have any names common with names of functional groups specified by 'FG'. Change the name of groups in either 'ID' or 'FG'")
  }
}

group_IDs <- function(data, prop, ID){
  # Species proportions to sum 
  props <- data[, prop] 
  
  # Unique IDs
  uIDs <- unique(ID) 
  
  # Number of unique ID terms
  nIDs <- length(uIDs) 
  
  sapply(uIDs, function(x){
    # Index of prop to sum
    idx <- which(ID == x) 
    if(length(idx) == 1){
      # No need to sum if only one element in group
      result <- props[, idx]
    } else {
      # Sum of species proportions
      result <-  rowSums(props[, idx])
    }
  })
}
#########################################

DI_data <- function(prop, FG, data, theta = 1, what = c("E","AV","FG","ADD","FULL")) {
  
  ########## RV change ###############
  if (inherits(data, 'tbl')){
    data <- as.data.frame(data)
  }
  
  ####################################
  
  if(length(what) == 1) {
    if("E" %in% what | "AV" %in% what) {
      if(length(prop) <= 2) {
        stop("You must have > 2 species to compute the 'E' or 'AV' variables") 
      }
    }
    if("ADD" %in% what) {
      if(length(prop) <= 3) {
        stop("You must have > 3 species to compute the 'ADD' variables") 
      }
    }
    result <- switch(EXPR = what,
                     "E" = DI_data_E_AV(prop = prop, data = data, theta = theta)$E,
                     "AV" = DI_data_E_AV(prop = prop, data = data, theta = theta)$AV,
                     "FG" = DI_data_FG(prop = prop, FG = FG, data = data, theta = theta)$FG,
                     "ADD" = as.matrix(DI_data_ADD(prop = prop, data = data, theta = theta)[[1]]),
                     "FULL" = DI_data_fullpairwise(prop = prop, data = data, theta = theta))
  } else {
    result <- list()
    if("E" %in% what) {
      if(length(prop) <= 2) {
        stop("You must have > 2 species to compute the 'E' or 'AV' variables") 
      }
      result$E <- DI_data_E_AV(prop = prop, data = data, theta = theta)$E
    }
    if("AV" %in% what) {
      if(length(prop) <= 2) {
        stop("You must have > 2 species to compute the 'E' or 'AV' variables") 
      }
      result$AV <- DI_data_E_AV(prop = prop, data = data, theta = theta)$AV
    }
    if("FG" %in% what) {
      if(missing(FG)) {
        stop("Please supply argument 'FG' to compute functional group variables")
      }
      result$FG <- DI_data_FG(prop = prop, FG = FG, data = data, theta = theta)$FG
    }
    if("ADD" %in% what) {
      if(length(prop) <= 3) {
        stop("The 'ADD' variables are only computed for > 3 species cases as the 'ADD' model is not informative for the 2 or 3 species case") 
      }
      result$ADD <- as.matrix(DI_data_ADD(prop = prop, data = data, theta = theta)[[1]])
    }
    if("FULL" %in% what) {
      result$FULL <- DI_data_fullpairwise(prop = prop, data = data, theta = theta)
    }
  }
  return(result)
}