test_that("functional redundancy works", {
  suppressWarnings(suppressMessages(library(dplyr)))
  
  # Test functional redundancy works for all DI models
  data("Switzerland")
  
  m_ID <- DI(y = "yield", prop = 4:7, data = Switzerland,
             DImodel = "ID") %>% suppressMessages()
  m_AV <- DI(y = "yield", prop = 4:7, data = Switzerland,
             DImodel = "AV") %>% suppressMessages()
  m_ADD <- DI(y = "yield", prop = 4:7, data = Switzerland,
              DImodel = "ADD") %>% suppressMessages()
  m_FG <- DI(y = "yield", prop = 4:7, data = Switzerland,
             DImodel = "FG", FG = c("G", "G", "L", "L"),
             estimate_theta = TRUE) %>% suppressMessages()
  m_FULL <- DI(y = "yield", prop = 4:7, data = Switzerland,
               DImodel = "FULL", theta = 1) %>% suppressMessages()
  
  m_treat_AV <- DI(y = "yield", prop = 4:7, data = Switzerland,
                   DImodel = "AV", treat = "nitrogen",
                   estimate_theta = TRUE) %>% suppressMessages()
  m_treat_ADD <- DI(y = "yield", prop = 4:7, data = Switzerland,
                    DImodel = "ADD", treat = "density", 
                    theta = 0.5) %>% suppressMessages()
  m_treat_FG <- DI(y = "yield", prop = 4:7, data = Switzerland,
                   DImodel = "FG", FG = c("G", "G", "L", "L"),
                   treat = "nitrogen", density = "density") %>%
    suppressMessages()
  
  # ID model
  expect_equal(fr_test_DI(m_ID, redundant = c("p1", "p2"),
                          delta = 2),
               fr_test(
                 fr_communities(prop = paste0("p", 1:4),
                                redundant = c("p1", "p2")),
                 object = m_ID, delta = 2
               ))
  
  # AV model
  AV_res <- fr_test_DI(m_AV, redundant = c("p3", "p4"),
                       delta = delta_sd(0.75))
  expect_equal(
    data.frame(AV_res)[, c("Set", "Estimate", "p_tost", "Equivalent")],
    data.frame(Set = c("Monos", "Within", "Between", "Between"),
               Estimate = c(3.484501, 3.567996, 1.742251, 1.742251),
               p_tost = c(0.9949681, 0.9999998, 0.8677154, 0.8677154),
               Equivalent = c(FALSE, FALSE, FALSE, FALSE)),
    tolerance = 1e-6, ignore_attr = TRUE
  )
  
  # FG model
  FG_res <- fr_test_DI(m_FG, redundant = c("p1", "p2"),
                       delta = delta_sd(0.75), alpha.level = 0.01)
  expect_equal(
    data.frame(FG_res)[, c("Set", "Estimate", "CI_level", "p_tost", "Equivalent")],
    data.frame(Set = c("Monos", "Within", "Between", "Between"),
               Estimate = c(-0.05667834, 1.80623352, -0.02833917, -0.02833917),
               CI_level = 0.98,
               p_tost = c(0.075838624, 0.750081932, 0.002077635, 0.002077635),
               Equivalent = c(FALSE, FALSE, TRUE, TRUE)),
    tolerance = 1e-6, ignore_attr = TRUE
  )
  
  # FULL model (delta as a number)
  FULL_res <- fr_test_DI(m_FULL, redundant = c("p1", "p2"),
                         delta = 5)
  expect_equal(
    data.frame(FULL_res)[, c("Set", "Estimate", "p_tost", "Equivalent")],
    data.frame(Set = c("Monos", "Within", "Between", "Between"),
               Estimate = c(-0.2493440, 2.5102927, 0.3163857, 0.1830209),
               p_tost = c(1.245899e-06, 7.672414e-03, 3.677897e-04, 2.672461e-04),
               Equivalent = c(TRUE)),
    tolerance = 1e-6, ignore_attr = TRUE
  )
  
  # ADD model (asymmetric interval)
  ADD_res <- fr_test_DI(m_ADD, redundant = c("p1", "p2"),
                        delta = c(-2, 3))
  expect_equal(
    data.frame(ADD_res)[, c("Set", "Estimate", "p_tost", "Equivalent")],
    data.frame(Set = c("Monos", "Within", "Between", "Between"),
               Estimate = c(-0.2493440, 4.6845570, 0.2497033, 0.2497033),
               p_tost = c(0.042154194, 0.968386952, 0.009417425, 0.009417425),
               Equivalent = c(TRUE, FALSE, TRUE, TRUE)),
    tolerance = 1e-6, ignore_attr = TRUE
  )
  
  # Models with treatment
  treat_res <- fr_test_DI(m_treat_ADD, redundant = c("p2", "p3"),
                          delta = delta_sd(1))
  expect_equal(
    data.frame(treat_res)[, c("Set", "Estimate", "p_tost", "Equivalent")],
    data.frame(Set = c("Monos", "Within", "Between", "Between"),
               Estimate = c(-6.724296, 2.660794, -3.628652, -3.628652),
               p_tost = c(0.9999773, 0.8583287, 0.9795225, 0.9795225),
               Equivalent = c(FALSE)),
    tolerance = 1e-6, ignore_attr = TRUE
  )
  # Check density was automatically added
  expect_equal(colnames(attr(treat_res, "contrast_matrix")),
               c("p1_ID", "p2_ID", "p3_ID", "p4_ID", 
                 "p1_add", "p2_add", "p3_add", "p4_add", "densitylow"))
  
  
  treat_res <- fr_test_DI(m_treat_FG, redundant = c("p1", "p2"),
                          delta = delta_sd(2))
  expect_equal(
    data.frame(treat_res)[, c("Set", "Estimate", "p_tost", "Equivalent")],
    data.frame(Set = c("Monos", "Within", "Between", "Between"),
               Estimate = c(-0.05667834, 2.37926953, -0.02833917, -0.02833917),
               p_tost = c(8.976175e-05, 2.558656e-01, 2.051899e-11, 2.051899e-11),
               Equivalent = c(TRUE, FALSE, TRUE, TRUE)),
    tolerance = 1e-6, ignore_attr = TRUE
  )
  # Check columns were added automatically
  expect_equal(colnames(attr(treat_res, "contrast_matrix")),
               c("p1_ID", "p2_ID", "p3_ID", "p4_ID",
                 "FG_bfg_G_L", "FG_wfg_G", "FG_wfg_L",
                 "nitrogen50", "densityhigh"))
  
  # Test if DI can be done with separate functions
  # Also methods work in a dplyr pipeline
  man_res <- fr_communities(paste0("p", 1:4), c("p3", "p4")) %>% 
    # Add columns
    dplyr::mutate(AV = DI_data(prop = paste0("p", 1:4), what = "AV", data =.),
           nitrogen = factor("150", c("50", "150"))) %>% 
    # fr_test
    fr_test(object = m_treat_AV, delta = 2, int_cols = "AV") %>% 
    suppressWarnings()
  
  expect_equal(
    data.frame(man_res)[, c("Set", "Estimate", "p_tost", "Equivalent")],
    data.frame(Set = c("Monos", "Within", "Between", "Between"),
               Estimate = c(3.484501, 1.862997, 1.742251, 1.742251),
               p_tost = c(0.9645218, 0.2598805, 0.2628668, 0.2628668),
               Equivalent = c(FALSE)),
    tolerance = 1e-6, ignore_attr = TRUE
  )
  expect_equal(attr(man_res, "contrast_matrix"),
               matrix(c(0, 0, 1.0, -1.0, 0.00, 0,
                        0, 0, 0.0,  0.0, 0.25, 0,
                        0, 0, 0.5, -0.5, 0.00, 0,
                        0, 0, 0.5, -0.5, 0.00, 0),
                      nrow = 4, byrow = TRUE), 
               ignore_attr = TRUE)
  
  # Ensure warning is returned if user specifies intearctions
  # with a DImodel
  expect_warning(fr_communities(paste0("p", 1:4), c("p3", "p4")) %>% 
                   # Add columns
                   mutate(AV = DI_data(prop = paste0("p", 1:4), what = "AV", data =.),
                          nitrogen = factor("150", c("50", "150"))) %>% 
                   # fr_test
                   fr_test(object = m_treat_AV, delta = 2, int_cols = "AV"),
                 "were already present in the")
  
  # Also test regular non-DI models work
  # With a base R example this time
  my_mod_treat <- lm(yield ~ 0 + p1 + p2 + p3 + p4 + int1 + int2 + int3 +
                       nitrogen + density,
                     data = Switzerland %>% 
                       dplyr::mutate(int1 = p1*p2, int2 = p2*p3, int3 = p2^p3))
  
  a <- fr_communities(paste0("p", 1:4), c("p1", "p2")) 
  
  a$int1 = a$p1*a$p2 
  a$int2 = a$p2*a$p3 
  a$int3 = a$p2^a$p3
  a$nitrogen = factor("150", levels = c("50", "150"))
  a$density = 100
  
  base_res <- fr_test(a, my_mod_treat, int_cols = 7:9,
                      delta = 1) 
  expect_equal(
    data.frame(base_res)$p_tost,
    c(0.4044036, 0.8503039, 0.8899725, 0.1196048),
    tolerance = 1e-6
  )

  ## Nine species example
  data("sim3")
  m9_FG <- DI(y = "response", prop = 4:12, 
              FG = c("G", "G", "G", "G", "G", "L", "L", "H", "H"),
              DImodel = "FG", theta = 0.25, data = sim3) %>% 
    suppressMessages() %>% suppressWarnings()
  
  nine_res_v1 <- fr_test_DI(m9_FG, redundant = c("p1", "p2", "p3", "p4", "p5"),
                            delta = 2)
  
  # Test all species
  nine_res_v2 <- fr_test_DI(m9_FG, redundant = paste0("p", 1:9),
                            delta = 2)
  
  # S3 methods work fine
  expect_snapshot(print(AV_res))
  expect_snapshot(summary(AV_res))
  # Because base R plot don't work with expect_snapshot
  expect_no_error(plot(AV_res))
  
  expect_snapshot(print(FULL_res))
  expect_snapshot(summary(FULL_res))
  # Because base R plot don't work with expect_snapshot
  expect_no_error(plot(FULL_res))
  
  expect_snapshot(print(base_res))
  expect_snapshot(summary(base_res))
  # Because base R plot don't work with expect_snapshot
  expect_no_error(plot(base_res))
  
  expect_snapshot(print(nine_res_v1))
  expect_snapshot(summary(nine_res_v1))
  # Because base R plot don't work with expect_snapshot
  expect_no_error(plot(nine_res_v1))
  
  expect_snapshot(print(nine_res_v2))
  expect_snapshot(summary(nine_res_v2))
  # Because base R plot don't work with expect_snapshot
  expect_no_error(plot(nine_res_v2))

})

test_that("functional redundancy throws appropriate errors", {
  
  suppressWarnings(suppressMessages(library(dplyr)))
  # Test functional redundancy works for all DI models
  data("Switzerland")
  
  m_AV <- DI(y = "yield", prop = 4:7, data = Switzerland,
             DImodel = "AV") %>% suppressMessages()
  m_simp <- lm(yield ~ 0 + p1 + p2 + p3 + p4,
               data = Switzerland)
  m_cust <- lm(yield ~ 0 + p1 + p2 + p3 + p4 +
                       nitrogen + density,
                     data = Switzerland)
  
  
  # Expected errors for fr_communities
  expect_error(fr_communities(), "`prop` should be")
  expect_error(fr_communities(prop = 1:4), "`redundant` should be")
  expect_error(fr_communities(prop = 1:4, redundant = 1), 
               "There should be at least two")
  expect_error(fr_communities(prop = 1:4, redundant = c(4, 5)), 
               "All names specified in `redundant` should be present in `prop`")
  
  # Expected errors for fr_contrasts
  expect_error(fr_contrasts(), "Please specify a regression model")
  expect_error(fr_contrasts(object = c("wrong")), "`data` is missing")
  expect_error(fr_contrasts(data = m_AV, object = Switzerland), 
               "`data` should be a data.frame returned by")
  expect_error(fr_contrasts(object = m_cust,
                       data = fr_communities(paste0("p", 1:4),
                                             c("p1", "p2"))),
               "`int_cols` is mandatory") %>% suppressWarnings()
  expect_error(fr_contrasts(data = fr_communities(paste0("p", 1:4),
                                                  c("p1", "p2")) %>% select(-Set),
                            object = m_simp, int_cols = NULL), 
               "The `Set` column in the data is important.")
  expect_error(fr_contrasts(data = fr_communities(paste0("p", 1:4),
                                                  c("p1", "p2")) %>% 
                              select(-Comm_ID),
                            object = m_simp, int_cols = NULL), 
               "The `Comm_ID` column in the data is important.")
  
  
  # Expected errors for fr_test_DI
  expect_error(fr_test_DI(), "Specify a fitted DI model")
  expect_error(fr_test_DI(c("wrong")), "Please provide a model fitted using the `DImodels`")
  expect_error(fr_test_DI(m_AV), "`redundant` should be")
  expect_warning(fr_test_DI(object = m_AV,
                            redundant = c("p1", "p2")),
                 "No `delta` supplied") 
  
  # Expected errors for fr_tost
  expect_error(fr_tost(), "Please specify a regression model")
  expect_error(fr_tost(c("wrong")), "Please specify a regression")
  expect_error(fr_tost(object = m_cust), "`data` is missing")
  expect_error(fr_tost(object = m_cust,
                       data = Switzerland),
               "`data` should be a data.frame returned by the")
  expect_warning(fr_tost(object = m_simp,
                         data = fr_contrasts(
                           fr_communities(paste0("p", 1:4),
                                          c("p1", "p2")),
                           m_simp, int_cols = NULL)
                         ),
                 "No `delta` supplied") 
  
  expect_error(fr_tost(object = m_simp,
                         data = fr_contrasts(
                           fr_communities(paste0("p", 1:4),
                                          c("p1", "p2")),
                           m_simp, int_cols = NULL)[, 1:3],
                         delta = 2),
               "`data` has 3 columns but object has 4 coefficients") 
  
  
  # Expected errors for fr_test
  expect_error(fr_test(), "Please specify a regression model")
  expect_error(fr_test(c("wrong")), "Please specify a regression")
  expect_error(fr_test(object = m_cust), "`data` is missing")
  expect_error(fr_test(object = m_cust,
                       data = fr_communities(paste0("p", 1:4),
                                             c("p1", "p2"))),
               "`int_cols` is mandatory") %>% suppressWarnings()
  expect_error(fr_test(object = m_cust,
                         data = fr_communities(paste0("p", 1:4),
                                               c("p1", "p2")),
                         int_cols = NULL,
                         delta = 2),
                 "The following variables required")
  expect_error(fr_test(object = m_cust,
                       data = fr_communities(paste0("p", 1:4),
                                             c("p1", "p2")) %>% 
                         mutate(density = "high",
                                nitrogen = "50"),
                       int_cols = NULL,
                       delta = 2),
               "The following variables were fitted as factors")
  expect_error(fr_test(object = m_cust,
                       data = fr_communities(paste0("p", 1:4),
                                             c("p1", "p2")),
                       int_cols = 8,
                       delta = 2),
               "Can't select columns past the end")
  
  expect_error(fr_test(object = m_cust,
                       data = fr_communities(paste0("p", 1:4),
                                             c("p1", "p2")),
                       int_cols = c("ola"),
                       delta = 2),
               "Can't select columns:")
  
  expect_error(fr_test(object = m_simp,
                       data = fr_communities(paste0("p", 1:4),
                                             c("p1", "p2")),
                       int_cols = NULL,
                       delta = -2),
               "If specifying a single value for `delta`, it should be a positive number.")
  expect_error(fr_test(object = m_simp,
                       data = fr_communities(paste0("p", 1:4),
                                             c("p1", "p2")),
                       int_cols = NULL,
                       delta = c(-3, -2)),
               "the first value should be negative while second should be positive")
  expect_error(fr_test(object = m_simp,
                       data = fr_communities(paste0("p", 1:4),
                                             c("p1", "p2")),
                       int_cols = NULL,
                       delta = c(3, -2)),
               "`delta_lo` must be lower than `delta_hi`")
  expect_error(fr_test(object = m_simp,
                       data = fr_communities(paste0("p", 1:4),
                                             c("p1", "p2")),
                       int_cols = NULL,
                       delta = c(1, 3, -2)),
               "Please specify `delta` as a positive finite vector of length one or two defining")
  
  
  expect_warning(fr_test(object = m_simp,
                         data = fr_communities(paste0("p", 1:4),
                                               c("p1", "p2")),
                         int_cols = NULL),
               "No `delta` supplied") 
  expect_no_error(fr_test(object = m_simp,
                          data = fr_communities(paste0("p", 1:4),
                                               c("p1", "p2")),
                          int_cols = NULL,
                          delta = 5))
  
  # Expected error in delta_sd
  expect_error(delta_sd(-5), "Please specify `k` as a single positive number.")
  
  # Expected errors in S3 methods
  expect_error(summary.fr_tost(m_simp),
               "`object` should be a data.frame")
  expect_error(print.fr_tost(m_simp),
               "`x` should be a data.frame")
  expect_error(plot.fr_tost(m_simp),
               "`x` should be a data.frame")
  expect_error(summary.fr_tost(fr_test(object = m_simp,
                                       data = fr_communities(paste0("p", 1:4),
                                                             c("p1", "p2")),
                                       int_cols = NULL,
                                       delta = 5),
                               max.print = 0), 
               "`max.print` must be a single positive number.")
})

