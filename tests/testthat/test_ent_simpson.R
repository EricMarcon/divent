# Combine all parameters
abundances <- paracou_6_abd[1, ]

testthat::test_that(
  "No estimator fails", {
    # Estimate diversity systematically
    ent_simpson.list <- lapply(
      # All estimators
      eval(formals(divent:::ent_simpson.numeric)$estimator), 
      function(estimator) {
        the_list <-lapply(
          # All probability estimators
          eval(formals(divent:::ent_simpson.numeric)$probability_estimator), 
          function(probability_estimator) {
            the_list <- lapply(
              # All unveilings
              eval(formals(divent:::ent_simpson.numeric)$unveiling), 
              function(unveiling) {
                # print(paste(estimator, probability_estimator, unveiling))
                suppressWarnings(
                  ent_simpson(
                    abundances,
                    estimator = estimator,
                    jack_alpha = 0.05,
                    jack_max = 10,
                    level = NULL,
                    probability_estimator = probability_estimator,
                    unveiling = unveiling,
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
    # Coerce to a dataframe
    ent_simpson.dataframe <- do.call(rbind, ent_simpson.list)
    
    # The min value must be UnveilJ / Chao2013 without unveiling
    testthat::expect_equal(
      min(ent_simpson.dataframe$entropy),
      ent_simpson(
        abundances, 
        estimator = "UnveilC",
        probability_estimator = "Chao2013",
        unveiling = "none"
      )$entropy
    )
  }
)


# Interpolation and extrapolation
sample_size <- abd_sum(abundances, as_numeric = TRUE)
levels <- c(sample_size / 2, round(sample_size * 1.5))

testthat::test_that(
  "No estimator fails during interpolation and extrapolation", {
    # Estimate diversity systematically
    ent_simpson.list <- lapply(
      # All probability estimators
      eval(formals(divent:::ent_simpson.numeric)$probability_estimator), 
      function(probability_estimator) {
        the_list <- lapply(
          # All unveilings
          eval(formals(divent:::ent_simpson.numeric)$unveiling), 
          function(unveiling) {
            the_list <-lapply(
              # All levels
              levels, 
              function(level) {
                # print(paste(estimator, probability_estimator, unveiling, level))
                suppressWarnings(
                  ent_simpson(
                    abundances,
                    estimator = "naive",
                    jack_alpha = 0.05,
                    jack_max = 10,
                    level = level,
                    probability_estimator = probability_estimator,
                    unveiling = unveiling,
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
    # Coerce to a dataframe
    ent_simpson.dataframe <- do.call(rbind, ent_simpson.list)
    
    # The min value must be at level 471, no matter the arguments
    testthat::expect_equal(
      min(ent_simpson.dataframe$entropy),
      as.numeric(
        unique(
          ent_simpson.dataframe[
            ent_simpson.dataframe$level == 471,
            "entropy"
          ]
        )
      )
    )
  }
)
