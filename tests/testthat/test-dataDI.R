## Test good scenarios
test_that("DI_data works", {
  data("sim1")
  
  ## Test E interaction structure
  exp_AV_output <- c(rep(0.24, 4), 0.375, rep(0.33, 6), rep(0, 4))
  names(exp_AV_output) <- seq(1, 60, 4)
  expect_equal(DI_data(prop = 3:6, data = sim1[seq(1, 60, 4), ], what = "AV"),
               exp_AV_output)
  
  ## Test AV interaction structure
  exp_E_output <- c(rep(0.64, 4), 1, rep(0.88, 6), rep(0, 4))
  names(exp_E_output) <- seq(1, 60, 4)
  expect_equal(DI_data(prop = 3:6, data = sim1[seq(1, 60, 4), ], what = "E"),
               exp_E_output)
  
  ## Test ADD interaction structure
  exp_ADD_output <- apply(sim1[seq(1, 60, 4), 3:6], 2, function(x) x*(1-x))
  colnames(exp_ADD_output) <- paste0("p", 1:4, "_add")
  expect_equal(DI_data(prop = 3:6, data = sim1[seq(1, 60, 4), ], what = "ADD"),
               exp_ADD_output)
  
  ## Test FG interaction structure
  exp_FG_output <- matrix(c(rep(0.16, 4), 0.25, 0.16, rep(0.25, 4), 0.16, rep(0, 4),
                            0.07, 0.07, 0.01, 0.01, 0.0625, 0.16, rep(0.04, 4), 0.01, rep(0, 4),
                            0.01, 0.01, 0.07, 0.07, 0.0625, 0.01, rep(0.04, 4), 0.16, rep(0, 4)),
                          ncol = 3)
  colnames(exp_FG_output) <- c("bfg_G_H", "wfg_G", "wfg_H")
  expect_equal(DI_data(prop = 3:6, data = sim1[seq(1, 60, 4), ], what = "FG",
                       FG = c("G", "G", "H", "H")),
               exp_FG_output)
  
  ## Test FULL interaction structure
  fmla <- as.formula("~ 0 + (p1 + p2 + p3 + p4)^2")
  exp_FULL_output <- model.matrix(fmla, data = sim1[seq(1, 60, 4), ])[, 5:10]
  expect_equal(DI_data(prop = 3:6, data = sim1[seq(1, 60, 4), ], what = "FULL"),
               exp_FULL_output)
  
  ## Test all interaction structures together
  exp_combined <- list("E" = exp_E_output,
                       "AV" = exp_AV_output,
                       "FG" = exp_FG_output,
                       "ADD" = exp_ADD_output,
                       "FULL" = exp_FULL_output)
  expect_equal(DI_data(prop = 3:6, data = sim1[seq(1, 60, 4), ],
                       FG = c("G", "G", "H", "H")),
               exp_combined)
  
})

## Test bad situations when function should throw error
test_that("DI_data gives appropriate error",{
  ## Ensure FG argument is passed when what includes "FG"
  data("sim1")
  expect_error(DI_data(prop = 3:6, data = sim1, what = "FG"))
  expect_error(DI_data(prop = 3:6, data = sim1, what = c("FULL", "FG")))
  
  ## Ensure proper names for Functional groups (reserved keywords can't be used)
  expect_error(DI_data(prop = 3:6, data = sim1, 
                       what = "FG", FG = c("g", "g", "in", "in")),
               regexp = "Please give your functional groups a different name")
  
  ## Ensure FG is specified as character strings
  expect_error(DI_data(prop = 3:6, data = sim1, 
                       what = "FG", FG = c(1, 1, 2, 2)),
               regexp = "FG argument takes character strings with functional group names")
  
  ## Ensure ADD is not calculated for data with less than four species
  data("sim0")
  expect_error(DI_data(prop = 3:5, data = sim0, what = "ADD"),
               regexp = "> 3 species")
  expect_error(DI_data(prop = 3:5, data = sim0, what = c("AV", "ADD")),
               regexp = "> 3 species")
  
  # AV and E can't be calculated for data with 2 species
  sim0$p4 <- 1 - sim0$p1
  expect_error(DI_data(prop = c("p1", "p4"), data = sim0, 
                       what = "AV"),
               regexp = "> 2 species")
  expect_error(DI_data(prop = c("p1", "p4"), data = sim0, 
                       what = c("AV", "FULL")),
               regexp = "> 2 species")
  expect_error(DI_data(prop = c("p1", "p4"), data = sim0, 
                       what = c("E", "FULL")),
               regexp = "> 2 species")
  
  # Data can't have column names ending with add
  dummy <- sim1
  colnames(dummy)[1:4] <- paste0("p", 1:4, "_add")
  expect_error(DI_data(prop = 3:6, data = dummy, what = "ADD"),
               regexp = "Certain column names cause internal conflicts when calculating additive interactions")
  
})

test_that("DI_data_prepare works", {
  data("Switzerland")
  # If y not provided a column of 0's will be added
  expect_warning(DI_data_prepare(data = Switzerland, prop = 4:7),
                 regexp = "y was not supplied, so a vector of zeros was created for the response")
  
  swiss_num <- Switzerland
  swiss_num$nitrogen <- as.numeric(swiss_num$nitrogen) 
  # If density provided as a numeric variable
  expect_warning(DI_data_prepare(data = swiss_num, prop = 4:7, y = 8,
                                 density = 2, treat = 3),
                 regexp = "not a factor")
  # If block specified as numeric variable
  expect_warning(DI_data_prepare(data = swiss_num, prop = 4:7, y = 8,
                                 block = 2, treat = 3),
                 regexp = "not a factor")
  
  # Don't accept NA values
  dummy <- rbind(Switzerland, c(NA, NA, NA, 0, 0, 0, 0, 0))
  expect_error(DI_data_prepare(data = dummy, prop = 4:7, y = 8,
                               treat = 3),
               regexp = "The dataset contains missing values. Please remove them prior to analysis")
  
  # Ensure species proportions sum to 1
  dummy <- rbind(Switzerland[, 4:8], c(0, 0.9, 0, 0, 0))
  expect_error(DI_data_prepare(data = dummy, prop = 1:4, y = 5),
               regexp = "One or more rows have species proportions that do not sum to 1. This must be corrected prior to analysis")
  
  # Ensure there are no negative proportions
  dummy <- rbind(Switzerland[, 4:8], c(-0.1, 0.9, 0.1, 0.1, 0))
  expect_error(DI_data_prepare(data = dummy, prop = 1:4, y = 5),
               regexp = "One or more rows have species proportions with values less than 0 or greater than 1. This must be corrected prior to analysis")
  
  
  # Warning if any species proportions sum to 0
  dummy <- rbind(Switzerland[, 4:8], c(0, 0, 0, 0, 0))
  expect_warning(DI_data_prepare(data = dummy, prop = 1:4, y = 5),
               regexp = "One or more rows in your dataset have ALL proportions equal to zero")
  
  
})
