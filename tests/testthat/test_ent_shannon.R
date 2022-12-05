# Combine all parameters
abundances <- paracou_6_abd[1, ]

testthat::test_that(
  "No estimator fails", {
    # Estimate diversity systematically
    ent_shannon.list <- lapply(
      # All estimators
      eval(formals(divent:::ent_shannon.numeric)$estimator), 
      function(estimator) {
        the_list <-lapply(
          # All probability estimators
          eval(formals(divent:::ent_shannon.numeric)$probability_estimator), 
          function(probability_estimator) {
            the_list <- lapply(
              # All unveilings
              eval(formals(divent:::ent_shannon.numeric)$unveiling), 
              function(unveiling) {
                # print(paste(estimator, probability_estimator, unveiling))
                suppressWarnings(
                  ent_shannon(
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
    ent_shannon.dataframe <- do.call(rbind, ent_shannon.list)
    
    # The min value must be UnveilJ / Chao2013 without unveiling
    testthat::expect_equal(
      min(ent_shannon.dataframe$entropy),
      ent_shannon(
        abundances, 
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
    ent_shannon.list <- lapply(
      # All probability estimators
      eval(formals(divent:::ent_shannon.numeric)$probability_estimator), 
      function(probability_estimator) {
        the_list <- lapply(
          # All unveilings
          eval(formals(divent:::ent_shannon.numeric)$unveiling), 
          function(unveiling) {
            the_list <-lapply(
              # All levels
              levels, 
              function(level) {
                # print(paste(probability_estimator, unveiling, level))
                suppressWarnings(
                  ent_shannon(
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
    ent_shannon.dataframe <- do.call(rbind, ent_shannon.list)
    
    # The min value must be UnveilJ / Chao2013 without unveiling
    testthat::expect_equal(
      min(ent_shannon.dataframe$entropy),
      ent_shannon(
        abundances, 
        probability_estimator = "Chao2013",
        unveiling = "none",
        level = 1413
      )$entropy
    )
  }
)
