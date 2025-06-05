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
  exp_pred <- c(13.2703349718639, 13.3043419781524, 17.4404149400718, 15.3497141686398)
  names(exp_pred) <- as.character(1:4)
  expect_equal(predict(mod, newdata = Switzerland[1:4, ]),
               exp_pred)
  
  # Ensure SE is accurate
  exp_se <- c(0.453997370504706, 0.453997370504706, 0.453997370504706, 0.453997370504706)
  names(exp_se) <- as.character(1:4)
  expect_equal(predict(mod, newdata = Switzerland[1:4, ], se.fit = TRUE)$se.fit,
               exp_se)
  
  # If only one row is specified
  expect_equal(suppressWarnings(predict(mod, newdata = Switzerland[4, ])),
               c("4" = 15.3497141686398))
  
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
  
  # If species don't sum to 1 then warning will be thrown
  expect_warning(predict(mod, newdata = data.frame(p1 = c(0.333, 1), p2 = c(0.333, 0),
                                                   p3 = c(0.333, 0), p4 = c(0, 1),
                                                   nitrogen = "50", density = "low")),
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
               c("1" = 12.8725915311115, "2" = 12.8424977926434, 
                 "3" = 17.0259686206005, "4" = 14.5222165084828))
  
  # Predict using type = "response"
  expect_equal(predict(mod_basic, newdata = Switzerland[1:4, ], type = "response"),
               c("1" = 12.8725915311115, "2" = 12.8424977926434, 
                 "3" = 17.0259686206005, "4" = 14.5222165084828))
  
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
                 regexp = "not supplied in newdata. Calculating the prediction")
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

test_that("Prediction works with grouped ID effects", {
  # Ensuring predictions work for all interaction structures
  data("Switzerland")
  swiss_num_treat <- Switzerland
  swiss_num_treat$nitrogen <- as.numeric(swiss_num_treat$nitrogen)
  
  # E model
  mod_int <- DI(y = "yield", prop = 4:7, treat = "nitrogen",
                density = "density", DImodel = "E",
                ID = c("I1","I1", "I1", "I1"),
                theta = 0.5,
                extra_formula = ~ p1:nitrogen, 
                data = swiss_num_treat)
  
  # If not all factor variables present in model are present in newdata warning will be thrown
  expect_equal(predict(mod_int, newdata = swiss_num_treat[1:3, ]),
               c(13.5089908, 15.7958394, 15.7958394),
               ignore_attr = TRUE)
  
  # AV model
  mod_int <- DI(y = "yield", prop = 4:7, treat = "nitrogen",
                density = "density", DImodel = "AV",
                ID = c("I1","I2", "I2", "I1"),
                estimate_theta = TRUE,
                extra_formula = ~ p1:nitrogen, 
                data = swiss_num_treat)
  
  # If not all factor variables present in model are present in newdata warning will be thrown
  expect_equal(predict(mod_int, newdata = swiss_num_treat[1:4, ]),
               c(13.2287655, 15.5557315, 15.5557315, 15.4019903),
               ignore_attr = TRUE)
  
  # ADD model
  mod_int <- DI(y = "yield", prop = 4:7, treat = "nitrogen",
                density = "density", DImodel = "ADD",
                ID = c("I1","I2", "I3", "I4"),
                estimate_theta = TRUE,
                extra_formula = ~ I1:nitrogen, 
                data = swiss_num_treat)
  
  # If not all factor variables present in model are present in newdata warning will be thrown
  expect_equal(predict(mod_int, newdata = swiss_num_treat[1:4, ]),
               c(13.485657, 13.449836, 17.641787, 15.021399),
               ignore_attr = TRUE)
  
  # FG model
  mod_int <- DI(y = "yield", prop = 4:7, treat = "nitrogen",
                density = "density", DImodel = "FG",
                FG = c("G", "G", "H", "H"),
                ID = c("I1","I1", "I2", "I2"),
                estimate_theta = TRUE,
                extra_formula = ~ p1:nitrogen, 
                data = Switzerland)
  
  # If not all factor variables present in model are present in newdata warning will be thrown
  expect_equal(predict(mod_int, newdata = Switzerland[1:6, ]),
               c(13.5117704, 13.5316355, 16.3710272, 16.3710272, 16.7358524, 14.2027248),
               ignore_attr = TRUE)
  
  # FULL model
  mod_int <- DI(y = "yield", prop = 4:7, treat = "nitrogen",
                density = "density", DImodel = "FULL",
                ID = c("I3","I3", "I3", "I3"),
                estimate_theta = TRUE,
                extra_formula = ~ p1:nitrogen, 
                data = swiss_num_treat)
  
  # If not all factor variables present in model are present in newdata warning will be thrown
  expect_equal(predict(mod_int, newdata = swiss_num_treat[1:4, ]),
               c(13.5538097, 14.3899227, 17.4391819, 15.1266995),
               ignore_attr = TRUE)
  
  # ID model
  mod_int <- DI(y = "yield", prop = 4:7, treat = "nitrogen",
                density = "density", DImodel = "ID",
                ID = c("I3","I3", "I3", "I3"),
                estimate_theta = TRUE, extra_formula = ~ p1:nitrogen, 
                data = swiss_num_treat)
  
  # If not all factor variables present in model are present in newdata warning will be thrown
  expect_equal(predict(mod_int, newdata = swiss_num_treat[1, ]),
               c(12.3940615),
               ignore_attr = TRUE)
})

# Test that predict function works with extra_formula and m
# missing variables
test_that("Prediction works with grouped ID effects", {
  # Ensuring predictions work for all interaction structures
  data("Switzerland")
  mod_FG <- DI(y = "yield", DImodel = "FG", FG = c("G", "G", "H", "H"),
               data = Switzerland, prop = 4:7,
               extra_formula = ~nitrogen*density)
  expect_equal(suppressWarnings(predict(mod_FG, newdata = Switzerland[1, 4:7])),
               c(11.939679),
               ignore_attr = TRUE)
  
  swiss_num_treat <- Switzerland
  swiss_num_treat$nitrogen <- as.numeric(swiss_num_treat$nitrogen)
  mod_FG <- DI(y = "yield", DImodel = "FG", FG = c("G", "G", "H", "H"),
               data = swiss_num_treat, prop = 4:7,
               extra_formula = ~nitrogen*density)
  expect_equal(suppressWarnings(predict(mod_FG, newdata = swiss_num_treat[1:2, 4:7])),
               c(13.036867, 13.070874),
               ignore_attr = TRUE)
})
  
# Testing CI and PI work

test_that("CI and PI work", {
  data("sim2")
  
  mod <- DI(y = "response", DImodel = "FULL",
            data = sim2, prop = 3:6)
  mod_lm <- lm(response ~ 0 + (p1 + p2 + p3 + p4)^2, data = sim2)
  
  # Base prediction in same
  expect_equal(predict(mod), predict(mod_lm))
  
  # CI is same
  expect_equal(predict(mod, interval = "conf"),
               predict(mod_lm, interval = "conf"))
  
  # PI is same
  expect_equal(predict(mod, interval = "pred", level = 0.9),
               suppressWarnings(predict(mod_lm, interval = "pred", level = 0.9)))
  
  
  # PI with se is same
  DI_pred <- predict(mod, interval = "pred", se.fit = TRUE)
  lm_pred <- suppressWarnings(predict(mod_lm, interval = "pred", se.fit = TRUE))
  expect_equal(as.numeric(DI_pred$se.fit), lm_pred$se.fit)
  expect_equal(DI_pred$fit, lm_pred$fit)
  
  # CI and PI with newdata work
  DI_pred <- predict(mod, newdata = sim2[1:5, ], interval = "conf", se.fit = TRUE)
  lm_pred <- suppressWarnings(predict(mod_lm, newdata = sim2[1:5, ], interval = "conf", se.fit = TRUE))
  expect_equal(as.numeric(DI_pred$se.fit), as.numeric(lm_pred$se.fit))
  expect_equal(DI_pred$fit, lm_pred$fit)
  
  DI_pred <- predict(mod, newdata = sim2[1:5, ], interval = "pred", se.fit = TRUE)
  lm_pred <- suppressWarnings(predict(mod_lm, newdata = sim2[1:5, ], interval = "pred", se.fit = TRUE))
  expect_equal(as.numeric(DI_pred$se.fit), as.numeric(lm_pred$se.fit))
  expect_equal(DI_pred$fit, lm_pred$fit)
  
  DI_pred <- predict(mod, newdata = sim2[1:5, ], interval = "none", se.fit = TRUE)
  lm_pred <- suppressWarnings(predict(mod_lm, newdata = sim2[1:5, ], interval = "none", se.fit = TRUE))
  expect_equal(as.numeric(DI_pred$se.fit), as.numeric(lm_pred$se.fit))
  expect_equal(DI_pred$fit, lm_pred$fit)
  
  # Ensure response type = "terms" works
  expect_equal(as.vector(predict(mod_lm, type = "terms")),
               as.vector(predict(mod, type = "terms")))
  
  DI_pred <- predict(mod, type = "terms", se.fit = TRUE)
  lm_pred <- predict(mod_lm, type = "terms", se.fit = TRUE)
  expect_equal(as.vector(DI_pred$fit), as.vector(lm_pred$fit))
  expect_equal(as.vector(DI_pred$se.fit), as.vector(lm_pred$se.fit))

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

# test fortify
# Test describe_model
test_that("fortify.DI works", {
  data(sim2)
  ## Fit model
  mod_FG <- DI(y = "response", FG = c("G", "G", "L", "L"), 
               prop = 3:6, data = sim2, DImodel = "FG")
  ## Describe model
  expect_equal(fortify(mod_FG),
               cbind(ggplot2:::fortify.lm(mod_FG), sim2[, 3:6])[, c(1:6, 13:16, 7:12)])
})
