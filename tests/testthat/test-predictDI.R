## Test good scenarios
test_that("predict function works", {
  data("Switzerland")
  
  ## Fit model
  mod <- DI(y = "yield", prop = 4:7, treat = "nitrogen",
            density = "density", DImodel = "AV",
            extra_formula = ~nitrogen:density,
            estimate_theta = TRUE, data = Switzerland)
  
  # If no newdata is specified
  expect_equal(predict(mod), mod$fitted.values)
  
  # If newdata is specified
  expect_equal(predict(mod, newdata = Switzerland[1:4, ]),
               c(13.2703349718639, 13.3043419781524, 17.4404149400718, 15.3497141686398))
  
  # Ensure SE is accurate
  expect_equal(predict(mod, newdata = Switzerland[1:4, ], se.fit = TRUE)$se.fit,
               c(0.457828639867846, 0.457828639867846, 0.457828639867846, 0.457828639867846))
  
  # If only one row is specified
  expect_equal(suppressWarnings(predict(mod, newdata = Switzerland[4, ])),
               c(15.3497141686398))
  
  # If not all factor variables present in model are present in newdata warning will be thrown
  expect_warning(predict(mod, newdata = Switzerland[1:4, -c(2)]),
                 regexp = "not supplied in newdata. Calculating for")
  
  # If any levels of factors in newdata were not present in model data error will be thrown
  expect_error(predict(mod, newdata = data.frame(p1 = c(0, 1), p2 = c(0, 0),
                                                 p3 = c(1, 0), p4 = c(0, 0),
                                                 nitrogen = c("150", "150"),
                                                 density = c("medium", "medium"))),
               regexp = "Predictions can't be made for these values.")
  
  # If not all species present in model warning will be thrown
  expect_warning(predict(mod, newdata = Switzerland[14:15, -c(4,5)]),
                 regexp = "were not present in newdata.")
  
  # If species don't sum to 1 then error will be thrown
  expect_error(predict(mod, newdata = Switzerland[12:13, -c(4,5)]),
                 regexp = "don't sum to 1")
  
  # If species approximately sum to 1, this will be allowed up to a tolerance level
  expect_error(predict(mod, newdata = data.frame(p1 = c(0.333, 1), p2 = c(0.333, 0),
                                                 p3 = c(0.333, 0), p4 = c(0, 0))),
               regexp = "don't sum to 1")
  
  
  # Prediction function works for custom formula
  mod_custom <- DI(custom_formula = "yield ~ p1 + p2 + p3 + p4",
                   data = Switzerland)
  expect_equal(suppressWarnings(predict(mod_custom, newdata = Switzerland[10:14, ])),
               c(13.1146310101112, 15.1826674910709, 11.0106935467002, 11.0673718905144, 17.9608268270468),
               ignore_attr = TRUE)
  
  # Predict for model without extra_formula
  mod_basic <- DI(y = "yield", prop = paste0("p", 1:4), 
                  DImodel = "FULL",
                  data = Switzerland)
  expect_equal(predict(mod_basic, newdata = Switzerland[1:4, ]),
               c(12.8725915311115, 12.8424977926434, 17.0259686206005, 14.5222165084828))
  
  # Predict using type = "response"
  expect_equal(predict(mod_basic, newdata = Switzerland[1:4, ], type = "response"),
               c(12.8725915311115, 12.8424977926434, 17.0259686206005, 14.5222165084828))
  
  # Predict for when additional treatment is numeric
  swiss_num_treat <- Switzerland
  swiss_num_treat$nitrogen <- as.numeric(swiss_num_treat$nitrogen)
  
  ## Fit model
  mod_num <- DI(y = "yield", prop = 4:7, treat = "nitrogen",
            density = "density", DImodel = "ADD",
            extra_formula = ~nitrogen:density, 
            data = swiss_num_treat)
  
  # If not numeric variables present in model are present in newdata warning will be thrown
  expect_warning(predict(mod_num, newdata = swiss_num_treat[1:4, -c(2)]),
                 regexp = "not supplied in newdata. Calculating for")
})

test_that("Prediction works for all interaction structures", {
  # Ensuring predictions work for all interaction structures
  data("Switzerland")
  swiss_num_treat <- Switzerland
  swiss_num_treat$nitrogen <- as.numeric(swiss_num_treat$nitrogen)
  
  # E model
  mod_int <- DI(y = "yield", prop = 4:7, treat = "nitrogen",
                density = "density", DImodel = "E",
                extra_formula = ~ p1:nitrogen, 
                data = swiss_num_treat)
  
  # If not all factor variables present in model are present in newdata warning will be thrown
  expect_equal(predict(mod_int, newdata = swiss_num_treat[1:4, ]),
               c(13.0437024719415, 12.9990123297592, 17.1350852916787, 15.0443845202466),
               ignore_attr = TRUE)
  
  # AV model
  mod_int <- DI(y = "yield", prop = 4:7, treat = "nitrogen",
                density = "density", DImodel = "AV",
                extra_formula = ~ p1:nitrogen, 
                data = swiss_num_treat)
  
  # If not all factor variables present in model are present in newdata warning will be thrown
  expect_equal(predict(mod_int, newdata = swiss_num_treat[1:4, ]),
               c(13.0437024719415, 12.9990123297592, 17.1350852916787, 15.0443845202466),
               ignore_attr = TRUE)
  
  # ADD model
  mod_int <- DI(y = "yield", prop = 4:7, treat = "nitrogen",
                density = "density", DImodel = "ADD",
                extra_formula = ~ p1:nitrogen, 
                data = swiss_num_treat)
  
  # If not all factor variables present in model are present in newdata warning will be thrown
  expect_equal(predict(mod_int, newdata = swiss_num_treat[1:4, ]),
               c(13.1059549340488, 13.0843467119269, 17.267817539884, 14.7640654277663),
               ignore_attr = TRUE)
  
  # FG model
  mod_int <- DI(y = "yield", prop = 4:7, treat = "nitrogen",
                density = "density", DImodel = "FG",
                FG = c("G", "G", "H", "H"),
                extra_formula = ~ p1:nitrogen, 
                data = swiss_num_treat)
  
  # If not all factor variables present in model are present in newdata warning will be thrown
  expect_equal(predict(mod_int, newdata = swiss_num_treat[1:4, ]),
               c(13.0895912595161, 13.1104305294664, 17.0564317980378, 14.9657310266058),
               ignore_attr = TRUE)
  
  # FULL model
  mod_int <- DI(y = "yield", prop = 4:7, treat = "nitrogen",
                density = "density", DImodel = "FULL",
                extra_formula = ~ p1:nitrogen, 
                data = swiss_num_treat)
  
  # If not all factor variables present in model are present in newdata warning will be thrown
  expect_equal(predict(mod_int, newdata = swiss_num_treat[1:4, ]),
               c(13.1059549340488, 13.0843467119269, 17.267817539884, 14.7640654277663),
               ignore_attr = TRUE)
  
  # ID model
  mod_int <- DI(y = "yield", prop = 4:7, treat = "nitrogen",
                density = "density", DImodel = "ID",
                extra_formula = ~ p1:nitrogen, 
                data = swiss_num_treat)
  
  # If not all factor variables present in model are present in newdata warning will be thrown
  expect_equal(predict(mod_int, newdata = swiss_num_treat[1:4, ]),
               c(12.5973542522238, 12.5526641100415, 16.688737071961, 14.598036300529),
               ignore_attr = TRUE)
})


test_that("contrasts function works", {
  data("Switzerland")
  
  ## Fit model
  mod <- DI(y = "yield", prop = 4:7, treat = "nitrogen",
            density = "density", DImodel = "AV",
            extra_formula = ~nitrogen:density,
            estimate_theta = TRUE, data = Switzerland)
  
  # Model should be a DImodels object
  expect_error(contrasts_DI(lm(yield ~ p1 + p2, data = Switzerland)), 
               regexp = "Please provied a DImodels model object")
  
  # Mandatory to specify either constrast_vars or contrast
  expect_error(contrasts_DI(mod), 
               regexp = "Provide either one of contrast_vars or constrast")
  
  # Can't specify both contrast_vars and contrast
  expect_error(contrasts_DI(mod, contrast_vars = 0, contrast = 0),
               regexp = "Provide only one of contrast_vars or constrast")
  
  # Ensure contrast vector is of appropriate type
  expect_error(contrasts_DI(mod, contrast = c("1", "-1", "0", "0")),
               regexp = "Specify contrast as either a numeric vector, list or matrix")
  
  # Ensure contrast vector is of proper length
  expect_error(contrasts_DI(mod, contrast = c(1, -1, 0, 0)),
               regexp = "Number of elements in contrasts vector should be a multiple of number of coefficients in model")
  
  # Ensure contrast has appropriate columns if specified as a matrix
  expect_error(contrasts_DI(mod, contrast = matrix(c(1, -1, 0, 0), ncol = 4)),
               regexp = "Number of columns in contrast matrix should be same as number of coefficients in model")
  
  # Ensure contrast has appropriate length if specified as list 
  expect_error(contrasts_DI(mod, contrast = list(1, -1, 0, 0)),
               regexp = "Lengths of each element of contrasts list should be same as number of coefficients in model")
  
  # Ensure contrast_vars are specified as a list
  expect_error(contrasts_DI(mod, contrast_vars = c("density" = c(-0.25, 0.25, 0.25, -0.25))),
               regexp = "Contrast variables should be specified as a nested list")
  
  # Ensure user specifies variables present in model in contrast_vars
  expect_error(contrasts_DI(mod, contrast_vars = list(p5 = c(0, 1))),
               regexp = "not present in model")
  
  # Ensure number of elements in contrast_vars are same as number of levels of factor in model
  expect_error(contrasts_DI(mod, contrast_vars = list("density" = c(-1, 1, 1))),
               regexp = "Lengths of each element of contrasts list should be same as levels of variable in model")
  
  # Correct examples
  the_C <- matrix(c(1, 1, -1, -1, 0, 0, 0, 0), nrow = 1)
  colnames(the_C) <- names(mod$coefficients[1:8])
  # Contrast as matrix
  expect_equal(contrasts_DI(mod, contrast = the_C),
               multcomp::glht(mod, linfct = the_C , 
                              coef = mod$coef[1:8], vcov = vcov(mod)),
               ignore_attr = TRUE)
  
  # Contrast as vector
  expect_equal(contrasts_DI(mod, contrast = c(1, 1, -1, -1, 0, 0, 0, 0)),
               multcomp::glht(mod, linfct = the_C , 
                              coef = mod$coef[1:8], vcov = vcov(mod)),
               ignore_attr = TRUE)
  
  # Contrast as list
  expect_equal(contrasts_DI(mod, contrast = list(c(1, 1, -1, -1, 0, 0, 0, 0))),
               multcomp::glht(mod, linfct = the_C , 
                              coef = mod$coef[1:8], vcov = vcov(mod)),
               ignore_attr = TRUE)
  
  # Using contrast_vars
  the_C <- matrix(c(1, 0, 0, 0, 0, 0, 0, 0), nrow = 1)
  colnames(the_C) <- names(mod$coefficients[1:8])
  expect_equal(contrasts_DI(mod, contrast_vars = list("p1" = c(1))),
               multcomp::glht(mod, linfct = the_C , 
                              coef = mod$coef[1:8], vcov = vcov(mod)),
               ignore_attr = TRUE)
  
})