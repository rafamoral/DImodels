## Test good scenarios
test_that("DI works with \"E\" model", {
  data("Switzerland")

  ## E model works
  mod_E <- DI(y = "yield", prop = 4:7, 
              density = "density", DImodel = "E",
              estimate_theta = TRUE, data = Switzerland)
  expect_equal(mod_E$coefficients,
               c(8.31251807959244, 8.36919642340669, 15.2626513599391, 11.7781500742191, 3.55873211379141, -0.137755987617646, 0.83571275568028),
               ignore_attr = TRUE)
  
  mod_E_theta <- DI(y = "yield", prop = 4:7, treat = "nitrogen",
              block = "density", DImodel = "E",
              estimate_theta = TRUE, data = Switzerland)
  expect_equal(mod_E_theta$coefficients,
               c(7.71767135989568, 7.77434970370992, 14.6678046402423, 11.1833033545223, 2.79449615220446, -0.137755987617647, 0.965375640775469, 0.760301181253625),
               ignore_attr = TRUE)
  
  mod_E_theta <- DI(y = "yield", prop = 4:7, treat = "nitrogen",
                     density = "density", DImodel = "E",
                     estimate_theta = TRUE, data = Switzerland)
  expect_equal(mod_E_theta$coef,
               c(8.5452910, 8.6019694, 15.4954243, 12.0109230, 2.794496, -0.9653756, 0.1377560, 0.7603012),
               ignore_attr = TRUE)

  mod_E_theta_custom <- DI(y = "yield", prop = 4:7, treat = "nitrogen",
                            block = "density", DImodel = "E",
                            theta = 0.5, data = Switzerland)
  expect_equal(mod_E_theta_custom$coef,
               c(7.71483343, 7.77151177, 14.66496671, 11.18046542, 1.276282, -0.13775599, 1.16877634, 0.500000),
               ignore_attr = TRUE)

  # Additional special conditions to cover all situations in code
  mod_E <- DI(y = "yield", prop = 4:7, treat = "nitrogen",
               density = "density", DImodel = "E",
               extra_formula = ~p1:nitrogen, data = Switzerland)
  expect_equal(mod_E$coef,
               c(8.5985566700200202, 8.52407309971615, 15.4175280362485392, 11.9330267505285228, 5.1554594166571261, -0.6406529039456219, 0.13775598761764632, -0.4003438816546424),
               ignore_attr = TRUE)

  mod_E <- DI(y = "yield", prop = 4:7, block = "nitrogen",
               density = "density", DImodel = "E",
               extra_formula = ~p1:nitrogen, data = Switzerland)
  expect_equal(mod_E$coef,
               c(8.5985566700200202, 8.52407309971615, 15.4175280362485392, 11.9330267505285228, 5.1554594166571261, -0.6406529039456219, 0.13775598761764632, -0.4003438816546424),
               ignore_attr = TRUE)

  mod_E <- DI(y = "yield", prop = 4:7, estimate_theta = T,
               DImodel = "E", extra_formula = ~p1:nitrogen, data = Switzerland)
  expect_equal(mod_E$coef,
               c(8.79816909765492, 8.32361781247475, 15.2170727490071, 11.7325714632871, 3.221081, -1.62146559961042, 0.804511561097006),
               ignore_attr = TRUE)

  mod_E <- DI(y = "yield", prop = 4:7, theta = 0.5,
               DImodel = "E", data = Switzerland)
  expect_equal(mod_E$coef,
               c(8.43218081507944, 8.48885915889368, 15.3823140954261, 11.8978128097061, 1.313556, 0.5),
               ignore_attr = TRUE)

})

# Most combinations are covered by the test cases in autoDI
# Checking for additional cases which weren't handled there
test_that("Additional test-cases for DI function", {
  data("Switzerland")
  available_models <- c("STR", "ID", "AV", "E",
                       "FG", "ADD", "FULL")
  
  # Model with treatment and density and extra formula
  coefficients <- list("STR" = c(14.7325931464327, 0.338889885147146, -0.238739909717057, -0.0357565936016199),
                       "ID" = c(12.7191573096262, 12.8599403616208, 19.8189815384772, 16.4000664930813, 0.0789205580693944, -0.368131033135327, -0.0480451360491652),
                       "AV" = c(8.37376681890422, 8.42410991369333, 15.3126245184171, 11.8231829008885, 13.8613027010167, -0.85707214528136, 0.17586229967076, 0.00361903522292134),
                       "E" = c(8.37376681890421, 8.42410991369333, 15.3126245184171, 11.8231829008885, 5.19798851288125, -0.857072145281361, 0.17586229967076, 0.00361903522292138),
                       "FG" = c(8.04275198699745, 8.09102486468054, 15.5573333454879, 12.066277338473, 18.2468914418737, 9.70477849123567, 0.697832803984883, -0.895087245635806, 0.188314586667893, 0.00480165465840326),
                       "ADD" = c(7.9211674672672, 8.16603686490083, 14.9115344207494, 12.6630338919457, 9.9818446117302, 8.4450637556713, 9.53731593345885, -0.127299299449124, -0.915704826082002, 0.195068114253214, 0.00544305113298685),
                       "FULL" = c(7.9211674672672, 8.16603686490083, 14.9115344207494, 12.6630338919457, 9.72985117122586, 24.0010539929582, 14.0697090606876, 22.1975434375367, 12.7996579039913, 0.712959437834094, -0.915704826082003, 0.195068114253215, 0.00544305113298689))
  models <- lapply(available_models, function(int_str){
    mod <- DI(y = "yield", prop = 4:7, 
              treat = "nitrogen", density = "density",
              extra_formula = ~ plot, FG = c("G", "G", "H", "H"),
              data = Switzerland, DImodel = int_str)
    expect_equal(mod$coefficients, coefficients[[int_str]], ignore_attr = TRUE)
  })
  
  # Model with treatment and extra formula
  coefficients <- list("STR" = c(14.3798837946354, 0.472226683397857, -0.0318349230648342),
                       "ID" = c(12.3938360952603, 12.5239439263854, 19.4746603897198, 16.0474206308018, -0.128420306719741, -0.0419468753200729),
                       "AV" = c(8.557347892295, 8.61333495504708, 15.5062508191916, 12.0212104610836, 13.7602667961902, -0.753432769542902, 0.000394896949268814),
                       "E" = c(8.55734789229499, 8.61333495504708, 15.5062508191915, 12.0212104610836, 5.16010004857132, -0.753432769542904, 0.000394896949268857),
                       "FG" = c(8.24122525582619, 8.29555168019982, 15.7631312061979, 12.276795854859, 18.1385234335331, 9.569598253824, 0.616277024715319, -0.783926834356985, 0.00134354297097426),
                       "ADD" = c(8.12851127157043, 8.37632785992546, 15.1254447318632, 12.8805635191488, 9.89883746064958, 8.38792815070444, 9.49008749067943, -0.164620580041107, -0.800463673245844, 0.00185799083219398),
                       "FULL" = c(8.12851127157042, 8.37632785992547, 15.1254447318632, 12.8805635191488, 9.58970841517837, 23.8708183990981, 13.949380629015, 22.0931793897904, 12.7052010184325, 0.62840971446269, -0.800463673245844, 0.00185799083219402))
  models <- lapply(available_models, function(int_str){
    mod <- DI(y = "yield", prop = 4:7, 
              treat = "nitrogen",
              extra_formula = ~ plot, FG = c("G", "G", "H", "H"),
              data = Switzerland, DImodel = int_str)
    expect_equal(mod$coefficients, coefficients[[int_str]], ignore_attr = TRUE)
  })

  # Model with just extra formula
  coefficients <- list("STR" = c(15.0069626872072, -0.0399466155307287),
                       "ID" = c(12.4310412685474, 12.5650242020063, 19.5187625341331, 16.0945446440076, -0.0441605423413487),
                       "AV" = c(8.8269711821793, 8.906469377341, 15.8177196107957, 12.351013621998, 13.5669928746415, -0.0130359273140788),
                       "E" = c(8.8269711821793, 8.906469377341, 15.8177196107957, 12.351013621998, 5.08762232799058, -0.0130359273140788),
                       "FG" = c(8.52915340159294, 8.60795970991156, 16.0806792043195, 12.6134336707359, 17.9368964812065, 9.25954556539232, 0.523075808493804, -0.0126406843101482),
                       "ADD" = c(8.42935873119389, 8.68891801506243, 15.4524555780312, 13.2219950563477, 9.68884597830234, 8.28101845613254, 9.42265158656329, -0.19258269370149, -0.0124262124953947),
                       "FULL" = c(8.42935873119389, 8.68891801506243, 15.4524555780312, 13.2219950563477, 9.27280723825924, 23.5933910126348, 13.7114270330073, 21.9188337911023, 12.5703292102002, 0.533011696686164, -0.0124262124953947))
  models <- lapply(available_models, function(int_str){
    mod <- DI(y = "yield", prop = 4:7, 
              extra_formula = ~ plot, FG = c("G", "G", "H", "H"),
              data = Switzerland, DImodel = int_str)
    expect_equal(mod$coefficients, coefficients[[int_str]], ignore_attr = TRUE)
  })
  
  # Model with density and extra formula
  coefficients <- list("STR" = c(15.2312041965283, -0.304786490592827, -0.0420291739082572),
                       "ID" = c(12.3324765061034, 12.4706928045875, 19.4277323846698, 16.0068157424997, 0.352692096189999, -0.0465788682945809),
                       "AV" = c(8.82706133670023, 8.90655213440679, 15.8177965992039, 12.3510848417485, 13.5672214100844, -0.000539836639348497, -0.013031701489559),
                       "E" = c(8.82706133670023, 8.90655213440679, 15.8177965992038, 12.3510848417485, 5.08770802878165, -0.00053983663934812, -0.013031701489559),
                       "FG" = c(8.52998003730568, 8.60871674423055, 16.0814216246409, 12.6141218147411, 17.939046732214, 9.26200409328415, 0.524917782616992, -0.0050756675955452, -0.0126009242479091),
                       "ADD" = c(8.43055794921086, 8.69006867899454, 15.4535466148142, 13.2230264659816, 9.6908851388531, 8.282631390706, 9.42410130360552, -0.191296194190481, -0.0075371760688523, -0.0123671497001649),
                       "FULL" = c(8.43055794921086, 8.69006867899454, 15.4535466148142, 13.2230264659816, 9.27645933338346, 23.5968798902277, 13.7147526930691, 21.921896442718, 12.5732286442846, 0.535747913239408, -0.0075371760688523, -0.0123671497001648))
  models <- lapply(available_models, function(int_str){
    mod <- DI(y = "yield", prop = 4:7, 
              density = "density", FG = c("G", "G", "H", "H"),
              extra_formula = ~ plot,
              data = Switzerland, DImodel = int_str)
    expect_equal(mod$coefficients, coefficients[[int_str]], ignore_attr = TRUE)
  })
})

## Special conditions for the ADD model along with theta_CI
test_that("ADD model works with theta", {
  data("Switzerland")
  
  mod_ADD_treat <- DI(y = "yield", treat = "nitrogen",
                      DImodel = "ADD", prop = 4:7, 
                      data = Switzerland,
                      estimate_theta = TRUE)
  expect_equal(mod_ADD_treat$coef,
              c(8.20728328146857, 8.44732385981388, 15.2043119670558, 13.0134282595101, 6.04716805325336, 5.09190694619538, 5.80285734859481, -0.936494203249058, -0.941661834080826, 0.78595066963828),
              ignore_attr = TRUE)
  #8.18542119397663, 8.4347651913239, 15.1857578074747, 12.9427523389733, 9.85581823286135, 8.35831709222467, 9.46561091111296, -0.183962680694284, -0.740738874359283
  mod_ADD <- DI(y = "yield", 
                DImodel = "ADD", prop = 4:7, 
                data = Switzerland,
                estimate_theta = TRUE)
  expect_equal(mod_ADD$coef,
              c(7.8513433527586, 8.09323716891178, 14.8497165412014, 12.6422322301889, 7.32032604290266, 6.22099797598873, 7.03400285269713, -0.545604749698444, 0.856025879261015),
              ignore_attr = TRUE)
  #7.88878936443855, 8.13813336178581, 14.8891259779366, 12.6461205094352, 10.1178646560904, 8.6203635154537, 9.72765733434199, 0.0780837425347466
  
  # Confidence interval for theta
  theta_CI(mod_ADD)
})

## Special conditions for the FG model
test_that("FG model works with theta", {
  data("Switzerland")
  
  mod_FG_treat <- DI(y = "yield", treat = "nitrogen",
                      DImodel = "FG", prop = 4:7,
                      FG = c("G", "G", "H", "H"),
                      data = Switzerland,
                      estimate_theta = TRUE)
  expect_equal(mod_FG_treat$coef,
               c(8.31620467951036, 8.3728830233246, 15.8716577347265, 12.3871564490065, 10.4737628188265, 3.88086134464133, -2.04718179494063, -0.972602380260069, 0.752389619705663),
               ignore_attr = TRUE)
  #8.18542119397663, 8.4347651913239, 15.1857578074747, 12.9427523389733, 9.85581823286135, 8.35831709222467, 9.46561091111296, -0.183962680694284, -0.740738874359283
  mod_FG <- DI(y = "yield", 
                DImodel = "FG", prop = 4:7, 
                FG = c("G", "G", "H", "H"),
                data = Switzerland,
                estimate_theta = TRUE)
  expect_equal(mod_FG$coef,
               c(7.94604929097243, 8.00272763478668, 15.4934048896756, 12.0089036039556, 12.8506780250522, 5.7074593892259, -1.04515673145053, 0.829932131271245),
               ignore_attr = TRUE)
  #7.88878936443855, 8.13813336178581, 14.8891259779366, 12.6461205094352, 10.1178646560904, 8.6203635154537, 9.72765733434199, 0.0780837425347466
  
})

## Ensure specific models can't be fit for models with less than four species
test_that("Correct number of species are present", {
  data("sim1")
  
  sim1$p5 <- 1 - sim1$p1
  sim1$p6 <- 1 - sim1$p1 - sim1$p2
  
  expect_error(DI(data = sim1, y = "response", 
                  prop = c("p1", "p5"), 
                  DImodel= "E"),
               "you must have > 2 species to fit model E")

  expect_error(DI(data = sim1, y = "response", 
                  prop = c("p1", "p5"), 
                  DImodel= "AV"),
               "you must have > 2 species to fit model AV")
  
  expect_error(DI(data = sim1, y = "response", 
                  prop = c("p1", "p5"), 
                  DImodel= "FG",
                  FG = c("G", "G")),
               "you must have > 2 species to fit model FG")
  
  expect_error(DI(data = sim1, y = "response", 
                  prop = c("p1", "p2", "p6"), 
                  DImodel= "ADD"),
               "you must have > 3 species to fit model ADD")

})

## Richness vs DI works
test_that("Richness vs DI works", {
  data("Switzerland")
  
  ## compare the richness model with DI alternatives
  t1 <- richness_vs_DI(y = "yield", prop = 4:7, data = Switzerland)
  m1 <- DI(y = "yield", prop = 4:7, data = Switzerland,
           DImodel = "AV", estimate_theta = TRUE)
  expect_equal(t1$coefficients, m1$coefficients)
    
  ## include the density effects in the linear predictors of the four models
  t2 <- richness_vs_DI(y = "yield", prop = 4:7, data = Switzerland, extra_formula = ~ density)
  m2 <- DI(y = "yield", prop = 4:7, data = Switzerland, DImodel = "AV",
           extra_formula = ~ density, estimate_theta = TRUE)
  expect_equal(t2$coefficients, m2$coefficients)
})

## S3 methods for class DI (extract, AIC, BIC)
test_that("S3 methods work", {
  data("Switzerland")
  
  mod <- DI(y = "yield", prop = 4:7, treat = "nitrogen", FG = c("G", "G", "H", "H"), 
            density = "density", DImodel = "FG", data = Switzerland)
  
  # Extract function
  int_terms <- DI_data(prop = 4:7, data = Switzerland, FG = c("G", "G", "H","H"))
  int_terms[3:5] <- lapply(int_terms[3:5], as.data.frame)
  expect_equal(extract(mod),
               int_terms,
               ignore_attr = TRUE)
  
  # AIC 
  expect_equal(AIC(mod), 260.57313407399)
  
  # AICc 
  expect_equal(AICc(mod), 264.43278)
  
  # BIC
  expect_equal(BIC(mod), 282.768211125751)
  
  # BICc
  expect_equal(BICc(mod), 290.911120732231)
  
})

## Custom formula
test_that("Custom formula example", {
  data("Switzerland")
  # Basic custom formula
  mod <- DI(custom_formula = "yield ~ p1 + p2 + p3 + p4",
            data = Switzerland)
  lm_ver <- stats::glm(yield ~ p1 + p2 + p3 + p4, 
                       data = Switzerland)
  expect_equal(mod$coef, 
               lm_ver$coef)
})

## Special scenarios which can throw errors
test_that("DI fails with appropriate message", {
  data("sim0")
  # To suppress rounding error warning
  sim0$p3 <- 1- sim0$p1 - sim0$p2
  
  # DImodel should be proper
  expect_error(DI(prop = 3:5, data = sim0, DImodel = "AV1"),
               regexp = "should be one of")
  
  # y can't be missing
  expect_error(DI(prop = 3:5, data = sim0),
               regexp = "You must supply a response variable name or column index through the argument 'y'")
  
  # Theta can't be negative
  expect_error(DI(y = "response", prop = 3:5, DImodel = "E",
                  data = sim0, theta = -1),
               regexp = "Please choose a positive value for theta")
  
  # Can't specify both custom and extra formula
  expect_error(DI(custom_formula = "response ~ p1 + p2 + p3",
                  extra_formula = "~ + p1:p2", data = sim0),
               regexp = "Please provide either custom_formula or extra_formula; not the two at the same time.")
  
  # Can't fit ADD model with 3 species
  expect_error(DI(y = "response", prop = 3:5, DImodel = "ADD",
                  data = sim0),
               regexp = "> 3 species")
  
  # Can't fit AV and E model with 2 species
  sim0$p4 <- 1 - sim0$p1
  expect_error(DI(y = "response", prop = c("p1", "p4"), 
                  DImodel = "AV", data = sim0),
               regexp = "> 2 species")
  expect_error(DI(y = "response", prop = c("p1", "p4"), 
                  DImodel = "E", data = sim0),
               regexp = "> 2 species")
  
  # FG should be specified when fitting FG model
  expect_error(DI(y = "response", prop = 3:5, DImodel = "FG",
                  data = sim0),
               regexp = "The argument FG must be specified alongside DImodel = 'FG'")
  
  data("Switzerland")
  # Can't estimate theta when using custom_formula
  expect_error(DI(custom_formula = "yield ~ p1 + p2 + p3 + p4",
                  data = Switzerland, estimate_theta = TRUE),
               regexp = "theta estimation not available when custom_formula is supplied")
  
  ## Warnings
  # Don't specify theta and estimate_theta, estimate_theta takes precedence 
  expect_warning(DI(y = "yield", prop = 4:7, DImodel = "AV",
                  data = Switzerland, theta = 0.5, estimate_theta = TRUE),
                regexp = "By specifying estimate_theta as TRUE, DI is overriding the specified theta value")
  
  # If custom_formula and DImodel both specified custom_formula takes precendence
  expect_warning(DI(y = "yield", prop = 4:7, DImodel = "FULL",
                    data = Switzerland, custom_formula = "yield ~ p1 + p2"),
                 regexp = "fitting custom DI model using supplied custom_formula instead of model")
  
})  

## Testing for reference model and model comparison along with anovaglm
test_that("Interior functions in DI_internal work", {
  data("Switzerland")
  
  mod <- DI(y = "yield", prop = 4:7, treat = "nitrogen", FG = c("G", "G", "H", "H"), 
            density = "density", DImodel = "FG", data = Switzerland,
            estimate_theta = TRUE)
  
  # DI_reference and DI_compare
  expect_equal(DI_reference,
               DI_compare(mod))
  
  # anovaDIglm
  obj <- anovaDIglm(mod)
  expect_equal(obj$`Resid. Dev`,
               c(13310.9015415132, 9103.40986879718, 6326.62608653603, 2402.17197250788, 407.144041330928, 145.443529829053, 140.640221756244, 139.663529974278, 127.40531115179, 127.082707045673),
               ignore_attr = TRUE)
  
  # Rao test
  obj <- anovaDIglm(mod, test = "Rao")
  expect_equal(obj$`Rao`,
               c(NA, 4207.491672716, 2776.78378226115, 3924.45411402815, 1995.02793117695, 261.700511501874, 4.80330807280922, 0.976691781966224, 12.258218822488, 0.322604106116728),
               ignore_attr = TRUE)
  
  ## anova_glmList
  mod1 <- DI(y = "yield", prop = 4:7, treat = "nitrogen",  
             density = "density", DImodel = "AV", 
             data = Switzerland,
             estimate_theta = TRUE)
  mod2 <- DI(y = "yield", prop = 4:7, treat = "nitrogen",  
             density = "density", DImodel = "FULL", 
             data = Switzerland,
             estimate_theta = TRUE)
  
  expect_equal(anova_glmlist(list(mod1, mod, mod2))$Deviance,
               c(NA, 37.8055122016968, 11.4973620431856))
  
  expect_equal(anova_glmlist(list(mod1, mod, mod2), test = "Rao")$Rao,
               c(NA, 37.8141225527972, 11.5427716510587))
})

# Test describe_model
test_that("describe_model works", {
  data(sim2)
  ## Fit model
  mod_FG <- DI(y = "response", FG = c("G", "G", "L", "L"), 
               prop = 3:6, data = sim2, DImodel = "FG")
  ## Describe model
  expect_equal(describe_model(mod_FG),
               "This model has 4 species with the Functional group interaction structure with FG values of (G, G, L, L).")
  
  
  mod_FULL <- DI(y = "response", estimate_theta = TRUE, 
                 prop = 3:6, data = sim2, DImodel = "FULL")
  expect_equal(describe_model(mod_FULL),
               "This model has 4 species with the Full pairwise interaction structure. The interaction terms have an associated theta (exponent on the interactions) value of '0.45'.")
  
  mod_AV <- DI(y = "response", ID = c("ID1", "ID1", "ID2", "ID2"),
               estimate_theta = TRUE, 
               prop = 3:6, data = sim2, DImodel = "AV")
  expect_equal(describe_model(mod_AV),
               "This model has 4 species with the Average interaction structure. The interaction terms have an associated theta (exponent on the interactions) value of '0.45'.The species identity effects are grouped as (ID1, ID1, ID2, ID2).")
  
  mod_STR <- DI(y = "response", FG = c("G", "G", "L", "L"), 
               prop = 3:6, data = sim2, DImodel = "STR")
  ## Describe model
  expect_equal(describe_model(mod_STR),
               "This model doesn't have any species identity or interaction effects and only consists of experimental structures.")
  
  mod_CUST <- DI(custom_formula = response ~ block, data = sim2)
  ## Describe model
  expect_equal(describe_model(mod_CUST),
               "This is a custom DI model.")
  
  # Proper error is thrown
  expect_error(describe_model(lm(1 ~ 1)),
               regexp = "`model` should be regression object of class <DI>")
})
