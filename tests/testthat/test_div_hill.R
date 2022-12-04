# Combine all parameters
abundances <- paracou_6_abd[1, ]
# integer and non-integer q's
orders <- (0:6) / 2

testthat::test_that(
  "No estimator fails", {
    # Estimate diversity systematically
    div_hill.list <- lapply(
      # All estimators
      eval(formals(divent:::div_hill.numeric)$estimator), 
      function(estimator) {
        the_list <-lapply(
          # All probability estimators
          eval(formals(divent:::div_hill.numeric)$probability_estimator), 
          function(probability_estimator) {
            the_list <-lapply(
              # All richness estimators
              eval(formals(divent:::div_hill.numeric)$richness_estimator), 
              function(richness_estimator) {
                the_list <-lapply(
                  # All q's
                  orders, 
                  function(q) {
                    the_list <- lapply(
                      # All unveilings
                      eval(formals(divent:::div_hill.numeric)$unveiling), 
                      function(unveiling) {
                        print(paste(estimator, probability_estimator, unveiling, richness_estimator, q))
                        suppressWarnings(
                          div_hill(
                            abundances,
                            q = q,
                            estimator = estimator,
                            level = NULL,
                            probability_estimator = probability_estimator,
                            unveiling = unveiling,
                            richness_estimator = richness_estimator,
                            as_numeric = FALSE,
                            check_arguments = TRUE
                          )
                        )
                      }
                    ) 
                    # Make a dataframe with the list to avoid nested lists
                    the_df <- do.call(rbind, the_list)
                  }
                )
                # Make a dataframe with the list to avoid nested lists
                the_df <- do.call(rbind, the_list)
              }
            )
            # Make a dataframe with the list to avoid nested lists
            the_df <- do.call(rbind, the_list)
          }
        )
        # Make a dataframe with the list to avoid nested lists
        the_df <- do.call(rbind, the_list)
      }
    )
    # Coerce to a dataframe
    div_hill.dataframe <- do.call(rbind, div_hill.list)
    
    # The min value must be the number of observed species
    testthat::expect_equal(
      min(div_hill.dataframe$diversity, na.rm = TRUE),
      div_hill(
        abundances, 
        q = max(orders), 
        probability_estimator = "Chao2013",
        unveiling = "none"
      )$diversity
    )
  }
)

# Interpolation and extrapolation
sample_size <- abd_sum(abundances, as_numeric = TRUE)
levels <- c(sample_size / 2, round(sample_size * 1.5))

testthat::test_that(
  "No estimator fails during interpolation and extrapolation", {
    # Estimate diversity systematically
    div_hill.list <- lapply(
      # All estimators
      eval(formals(divent:::div_hill.numeric)$estimator), 
      function(estimator) {
        the_list <-lapply(
          # All probability estimators
          eval(formals(divent:::div_hill.numeric)$probability_estimator), 
          function(probability_estimator) {
            the_list <-lapply(
              # All richness estimators
              eval(formals(divent:::div_hill.numeric)$richness_estimator), 
              function(richness_estimator) {
                the_list <-lapply(
                  # All levels
                  levels, 
                  function(level) {
                    the_list <- lapply(
                      # All q's
                      orders, 
                      function(q) {
                        the_list <- lapply(
                          # All unveilings
                          eval(formals(divent:::div_hill.numeric)$unveiling), 
                          function(unveiling) {
                            print(paste(estimator, probability_estimator, unveiling, richness_estimator, q, level))
                            suppressWarnings(
                              div_hill(
                                abundances,
                                q = q,
                                estimator = estimator,
                                level = NULL,
                                probability_estimator = probability_estimator,
                                unveiling = unveiling,
                                richness_estimator = richness_estimator,
                                as_numeric = FALSE,
                                check_arguments = TRUE
                              )
                            )
                          }
                        ) 
                        # Make a dataframe with the list to avoid nested lists
                        the_df <- do.call(rbind, the_list)
                      }
                    )
                    # Make a dataframe with the list to avoid nested lists
                    the_df <- do.call(rbind, the_list)
                  }
                )
                # Make a dataframe with the list to avoid nested lists
                the_df <- do.call(rbind, the_list)
              }
            )
            # Make a dataframe with the list to avoid nested lists
            the_df <- do.call(rbind, the_list)
          }
        )
        # Make a dataframe with the list to avoid nested lists
        the_df <- do.call(rbind, the_list)
      }
    )
    # Coerce to a dataframe
    div_hill.dataframe <- do.call(rbind, div_hill.list)
    
    # The min value must be the number of observed species
    testthat::expect_equal(
      min(div_hill.dataframe$diversity, na.rm = TRUE),
      div_hill(
        abundances, 
        q = max(orders), 
        probability_estimator = "Chao2013",
        unveiling = "none"
      )$diversity
    )
  }
)
