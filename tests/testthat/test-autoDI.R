## Test good scenarios
test_that("autoDI works", {
  # Model with theta significant (no experimental structures)
  data("sim2")
  mod <- autoDI(y = "response", prop = 3:6, data = sim2,
                selection = "AIC",
                FG = c("G", "G", "H", "H"), step0 = TRUE)
  expect_equal(mod$coef, 
               c(10.8661494474281, 9.85152282405146, 8.4693312656099, 8.1479838630125, 6.99549861541515, 0.453346431366962),
               ignore_attr = TRUE)
  
  # Model with theta significant and treatment
  mod <- autoDI(y = "response", prop = 3:6, data = sim2,
                FG = c("G", "G", "H", "H"), selection = "AIC",
                treat = "block", step0 = TRUE)
  expect_equal(mod$coef, 
               c(9.38374994183338, 8.36912331845678, 6.98693176001522, 6.66558435741782, 6.99552456900258, 1.49453333333333, 2.2542, 2.18086666666667, 0.453347686849901),
               ignore_attr = TRUE)
  
  
  # Model where theta not significant (treat and density)
  data("Switzerland")
  mod <- autoDI(y = "yield", prop = 4:7, treat = "nitrogen", 
     density = "density", data = Switzerland, step0 = TRUE)
  expect_equal(mod$coef, 
               c(7.95766735824737, 8.20701135559464, 14.9580039717454, 12.714998503244, 10.0411709753684, 24.3274154382015, 14.4111121470316, 22.5631845982022, 13.1803407057576, 1.1086838807011, -0.137755987617647),
               ignore_attr = TRUE)
  
  # Additional models with different combinations of treat, density and block to cover different possible conditions
  data("Switzerland")
  # Treat and block
  mod <- autoDI(y = "yield", prop = 4:7, treat = "nitrogen", 
                block = "density", data = Switzerland, step0 = TRUE,
                FG = c("G", "G", "H", "H"),
                selection = "BICc")
  expect_equal(mod$coef, 
               c(8.05400018501388, 8.11067852882813, 15.5787518803547, 12.0942505946347, 18.6205132222982, 10.0411709753684, 1.1086838807011, -0.137755987617647),
               ignore_attr = TRUE)
  
  # Density and block
  mod <- autoDI(y = "yield", prop = 4:7, density = "nitrogen", 
                block = "density", data = Switzerland, step0 = TRUE,
                FG = c("G", "G", "H", "H"),
                selection = "AICc")
  expect_equal(mod$coef, 
               c(7.41248029170225, 7.46915863551648, 14.9679333472748, 11.4834320615547, 10.4737627252811, 3.88086127696886, -2.04718182436897, -0.137755987617647, 0.972602383743152, 0.752389615880585),
               ignore_attr = TRUE)
  
  # Only treat
  mod <- autoDI(y = "yield", prop = 4:7, treat = "nitrogen", 
                data = Switzerland, step0 = TRUE,
                FG = c("G", "G", "H", "H"),
                selection = "BIC")
  expect_equal(mod$coef, 
               c(7.98512219120506, 8.0418005350193, 15.5098738865459, 12.0253726008259, 18.6205132222982, 10.0411709753684, 1.1086838807011),
               ignore_attr = TRUE)
  
  # Only density
  mod <- autoDI(y = "yield", prop = 4:7, density = "nitrogen", 
                data = Switzerland, step0 = TRUE,
                FG = c("G", "G", "H", "H"),
                selection = "BIC")
  expect_equal(mod$coef, 
               c(8.28175402074315, 8.33843236455739, 15.806505716084, 12.322004430364, 18.0964203758402, 9.51707812891037, 0.58459103424304, -0.740738874359283),
               ignore_attr = TRUE)
  
  # Only block
  mod <- autoDI(y = "yield", prop = 4:7,
                block = "density", data = Switzerland, step0 = TRUE,
                FG = c("G", "G", "H", "H"),
                selection = "BIC")
  expect_equal(mod$coef, 
               c(8.05400018501388, 8.11067852882813, 15.5787518803547, 12.0942505946347, 18.6205132222982, 10.0411709753684, 1.1086838807011, -0.137755987617647),
               ignore_attr = TRUE)
  
  # All three
  block_data <- rbind(Switzerland, Switzerland, Switzerland)
  block_data$block <- rep(1:68, times = 3)
  mod <- autoDI(y = "yield", prop = 4:7, treat = "nitrogen",
                density = "density", block = "block",
                data = block_data, step0 = TRUE,
                FG = c("G", "G", "H", "H"),
                selection = "BIC")
  expect_equal(mod$coef, 
               c(7.33967637998198, 7.56306535822724, 14.2972367312356, 12.0983242028807, 0.0223398871371149, -1.71985109738412, 0.372981858061387, 3.99219198545826, 13.9846441926885, 7.5898745228345, 12.786357226827, 6.86155086894384, -1.92535067925968, 0.727780561011365),
               ignore_attr = TRUE)
  
})

## Test special scenarios
test_that("Special conditions", {
  data("sim0")
  
  # Ensure lack of fit not possible when communities are not repeated
  model_data <- sim0[1:16, ]
  mod <- autoDI(y = "response", prop = 3:5, 
                data = model_data, selection = "AIC")
  expect_equal(mod$coef,
               c(24.9971103842124, 17.7232005181953, 16.0731754351011, 44.5461945024041, 14.2889751130691, 52.5908021532639),
               ignore_attr = TRUE)
  
  # Ensure response is specified
  expect_error(autoDI(), 
               regexp = "You must supply a response variable name or column index through the argument 'y'")
  
  # Ensure specific models are not fit with 2 designs having two species
  sim0$p4 <- 1 - sim0$p1
  mod <- autoDI(y = "response", prop = c("p1", "p4"), 
                data = sim0, selection = "AIC")
  expect_equal(mod$coef,
               c(24.1179563413357, 22.1033530065088, 24.6768038740904),
               ignore_attr = TRUE)
  
  # Drop ADD model if selection F-test and FG is specified
  expect_message(autoDI(FG = c("G","G","H"), data = sim0, y = "response", prop = 3:5),
                 regexp = "Since 'Ftest' was specified as selection criterion and functional groups were specified, dropping the ADD model as it is not nested within the FG model")
})
