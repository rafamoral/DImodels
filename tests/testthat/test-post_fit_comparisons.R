# compare_communities test
test_that("compare_communities works",{
  data(sim2)
  sim2$block <- as.factor(sim2$block)
  ## Fit model
  mod_ID <- DI(y = "response", FG = c("G", "G", "L", "L"), 
               prop = 3:6, estimate_theta = TRUE,
               data = sim2, DImodel = "ID", ID = c("G1", "G2", "L1", "L2"))
  mod_E <- DI(y = "response", FG = c("G", "G", "L", "L"), 
               prop = 3:6, data = sim2, DImodel = "E", theta = 0.7)
  mod_AV <- DI(y = "response", FG = c("G", "G", "L", "L"), 
               prop = 3:6, data = sim2, DImodel = "AV", theta = 0.7)
  mod_FG <- DI(y = "response", FG = c("G", "G", "L", "L"), 
               prop = 3:6, data = sim2, DImodel = "FG", estimate_theta = TRUE)
  mod_ADD <- DI(y = "response", FG = c("G", "G", "L", "L"), 
                prop = 3:6, data = sim2, DImodel = "ADD",
                ID = c("p1", "p2", "p3", "p4"), treat = "block")
  mod_FULL <- DI(y = "response", FG = c("G", "G", "L", "L"), 
                 prop = 3:6, data = sim2, DImodel = "FULL", extra_formula = ~block,
                 estimate_theta = TRUE, ID = c("p1", "p1", "p2", "L"))
  
  # Test cases
  test_data <- sim2[sim2$block == 1, 3:6]
  test_data <- test_data[c(1:5, 12:15), ]
  rownames(test_data) <- c("p1_dom", "p2_dom", "p3_dom", "p4_dom", "cent",
                           "p1", "p2", "p3", "p4")
  expect_equal(unname(compare_communities(mod_ID, test_data[-5, ])$CLD),
               rep("a", 8))
  expect_equal(unname(compare_communities(mod_AV, test_data[c(5:9), ], adjust = "none")$CLD),
               c("a", "b", "bc", "cd", "d"))
  expect_equal(unname(compare_communities(mod_E, test_data[c(5:9), ], adjust = "none")$CLD),
               c("a", "b", "bc", "cd", "d"))
  expect_equal(unname(compare_communities(mod_FG, test_data[c(5:9), ], adjust = "free")$CLD), 
               c("a", "b", "bc", "bc", "c"))
  expect_equal(unname(compare_communities(mod_ADD, `$<-`(test_data[1:5,], "block", "1"))$CLD), 
               c("a", "ab", "b", "b", "c"))
  expect_equal(unname(compare_communities(mod_FULL,`$<-`(test_data[1:5,], "block", "1"), adjust = "none")$CLD), 
               c("a", "b", "bc", "c", "d"))
  expect_equal(unname(compare_communities(mod_AV, `$<-`(test_data[5:9,], "AV", 0), ref = "cent", adjust = "none")$CLD), 
               c("a", "b", "a", "a", "b"))
  expect_equal(unname(compare_communities(mod_AV, test_data[c(1:5), ], ref = "cent", adjust = "none")$CLD), 
               c("a", "b", "b", "b", "b"))
  expect_equal(unname(compare_communities(mod_AV, test_data[c(1:5), ], ref = 5, adjust = "none")$CLD), 
               c("a", "b", "b", "b", "b"))
  # warnings
  expect_warning(compare_communities(mod_ADD, test_data[c(1:5), ]),
                 "not supplied in newdata.")
  sim2_block_num <- sim2
  sim2_block_num$block <- as.numeric(sim2_block_num$block)
  mod_ADD_num <- DI(y = "response", FG = c("G", "G", "L", "L"), 
                prop = 3:6, data = sim2_block_num, DImodel = "ADD",
                ID = c("p1", "p2", "p3", "p4"), treat = "block")
  expect_warning(compare_communities(mod_ADD_num, test_data[c(1:5), 1:4]),
                 "not supplied in newdata.")
  
  expect_warning(compare_communities(mod_AV, test_data[c(1:5), 2:4]),
                 "The species proportions for some communities do not sum to 1")
  expect_warning(contrast_matrix(mod_AV, list(test_data[c(1:5), ])),
                 "`contrast_vars` should be specified as a <data.frame> or <matrix> ")
  
  # errors
  expect_error(compare_communities(mod_AV, test_data[1,]),
               "The `data` should contain more than 1 community for comparisons.")
  expect_error(compare_communities(mod_AV, list(p1 = 1)),
               "Please provide a data.frame containing the species communities to be compared in `data`")
  mod_lm <- lm(response ~ p1 + p2 + p3, data = sim2)
  expect_error(compare_communities(mod_lm, newdata = test_data),
               "Please provide a DImodels model object in `object`.")
  
  expect_error(compare_communities(mod_AV, test_data, ref = "5"),
               "`ref` should be one of the values from the row.names of `data`.")
  expect_error(compare_communities(mod_AV, test_data, ref = 10),
               "If specified as a number, `ref` can take values between 1 and ")
  expect_error(compare_communities(mod_AV, test_data, ref = test_data),
               "ref must be a row name or row index")
  
  test_data2 <- test_data[1:5, ]
  test_data2$block <- factor(5, levels = 1:5)
  expect_error(compare_communities(mod_FULL, test_data2),
               "Predictions can't be made for these values.")
})

# Get row differences works
test_that("tests for row_difference", {
  test_data <- sim2[sim2$block == 1, 3:6]
  test_data <- test_data[c(1:5, 12:15), ]
  
  expect_error(get_row_differences(as.character(test_data)),
               "`data` should be a data.frame or a matrix")
  expect_error(get_row_differences(`$<-`(test_data, "p1", "5")),
               "All values in `data` should be numeric.")
  
  expect_equal(unname(unlist(get_row_differences(test_data[1:2, ]))),
               c(0.6,-0.6,0,0))
})

# Testing contrast function
test_that("contrasts function works", {
  data("Switzerland")
  
  ## Fit model
  mod <- DI(y = "yield", prop = 4:7, treat = "nitrogen",
            density = "density", DImodel = "AV",
            extra_formula = ~nitrogen:density,
            estimate_theta = TRUE, data = Switzerland)
  
  # Model should be a DImodels object
  expect_error(contrasts_DI(lm(yield ~ p1 + p2, data = Switzerland)), 
               regexp = "Please provide a DImodels model object")
  
  # Mandatory to specify either constrast_vars or contrast
  expect_error(contrasts_DI(mod), 
               regexp = "Provide either one of `contrast_vars` or `constrast`")
  
  # Can't specify both contrast_vars and contrast
  expect_warning(contrasts_DI(mod, contrast_vars = 0, contrast = matrix(0, ncol = 8)),
                 regexp = "Provide only one of `contrast_vars` or `constrast`")
  
  # Ensure contrast vector throws error if it's not a matrix or data-frame
  expect_error(contrasts_DI(mod, contrast = c( "1", "-1", "0", "0")),
               regexp = "Specify contrast as either a")
  
  # Ensure contrast has appropriate columns if specified as a matrix
  expect_error(contrasts_DI(mod, contrast = matrix(c(1, -1, 0, 0), ncol = 4)),
               regexp = "Number of columns in contrast matrix should be same as number of coefficients in model")
  
  # Ensure contrast has appropriate length if specified as list 
  expect_error(contrasts_DI(mod, contrast = list(1, -1, 0, 0)),
               regexp = "Lengths of each element of contrasts list should be same as number of coefficients in model")
  
  # Correct examples
  the_C <- matrix(c(1, 1, -1, -1, 0, 0, 0, 0), nrow = 1)
  colnames(the_C) <- names(mod$coefficients[1:8])
  # Contrast as matrix
  expect_equal(contrasts_DI(mod, contrast = the_C),
               multcomp::glht(mod, linfct = the_C , 
                              coef = mod$coef[1:8], vcov = vcov(mod)),
               ignore_attr = TRUE)
  
  # Contrast as vector
  expect_equal(contrasts_DI(mod, contrast = matrix(c(1, 1, -1, -1, 0, 0, 0, 0), nrow = 1)),
               multcomp::glht(mod, linfct = the_C , 
                              coef = mod$coef[1:8], vcov = vcov(mod)),
               ignore_attr = TRUE)
  
  # Contrast as list
  expect_equal(contrasts_DI(mod, contrast = list(c(1, 1, -1, -1, 0, 0, 0, 0))),
               multcomp::glht(mod, linfct = the_C , 
                              coef = mod$coef[1:8], vcov = vcov(mod)),
               ignore_attr = TRUE)
  
  # Using contrast_vars
  
  
})


# test add_Int_ID works
test_that("add_Int_ID works",{
  data(sim2)
  
  ## Fit model
  mod_ID <- DI(y = "response", FG = c("G", "G", "L", "L"), 
               prop = 3:6, estimate_theta = TRUE,
               data = sim2, DImodel = "ID", ID = c("G", "G", "L", "L"))
  mod_AV <- DI(y = "response", FG = c("G", "G", "L", "L"), 
               prop = 3:6, data = sim2, DImodel = "AV")
  mod_FG <- DI(y = "response", FG = c("G", "G", "L", "L"), 
               prop = 3:6, data = sim2, DImodel = "FG")
  mod_FG_treat <- DI(y = "response", FG = c("G", "G", "L", "L"), 
               prop = 3:6, data = sim2, treat = "block", DImodel = "FG")
  mod_ADD <- DI(y = "response", FG = c("G", "G", "L", "L"), 
                prop = 3:6, data = sim2, DImodel = "ADD",
                ID = c("p1", "p2", "p3", "p4"))
  mod_FULL <- DI(y = "response", FG = c("G", "G", "L", "L"), 
                 prop = 3:6, data = sim2, DImodel = "FULL",
                 estimate_theta = TRUE, ID = c("p1", "p1", "p2", "L"))
  
  # See all interaction and ID columns are added
  new_dat <- sim2[1:4, 3:6]
  expect_equal(colnames(add_int_ID(mod_ID, new_dat)),
               c(colnames(new_dat), "G", "L"))
  expect_equal(colnames(add_int_ID(mod_AV, new_dat)),
               c(colnames(new_dat), "AV", paste0("p", 1:4, "_ID")))
  expect_equal(colnames(add_int_ID(mod_FG, new_dat)),
               c(colnames(new_dat), "FG_", paste0("p", 1:4, "_ID")))
  expect_equal(colnames(add_int_ID(mod_FG_treat, new_dat)),
               c(colnames(new_dat), "FG_bfg_G_L", "FG_wfg_G", "FG_wfg_L", paste0("p", 1:4, "_ID")))
  expect_equal(colnames(add_int_ID(mod_ADD, new_dat)),
               c(colnames(new_dat), paste0("p", 1:4, "_add")))
  expect_equal(colnames(add_int_ID(mod_FULL, new_dat[1, ])),
               c(colnames(new_dat), "p1:p2","p1:p3","p1:p4","p2:p3","p2:p4","p3:p4",  "L"))
})
