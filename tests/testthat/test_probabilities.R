# Combine all parameters
abundances <- paracou_6_abd[1, ]

testthat::test_that(
  "No estimator fails", {
    # Estimate diversity systematically
    probabilities.list <- lapply(
      # All estimators
      eval(formals(divent:::probabilities.numeric)$estimator), 
      function(estimator) {
        the_list <-lapply(
          # All unveilings
          eval(formals(divent:::probabilities.numeric)$unveiling), 
          function(unveiling) {
          the_list <- lapply(
            # All richness estimators
            eval(formals(divent:::probabilities.numeric)$richness_estimator), 
            function(richness_estimator) {
              the_list <- lapply(
                # All coverage estimators
                eval(formals(divent:::probabilities.numeric)$coverage_estimator), 
                function(coverage_estimator) {
                  # print(paste(estimator, unveiling, richness_estimator, coverage_estimator))
                  # Forbidden combination raises an error
                  if ((richness_estimator == "rarefy" & unveiling == "none")) {
                    NULL
                  } else {
                    suppressWarnings(
                      probabilities(
                        abundances,
                        estimator = estimator,
                        unveiling = unveiling,
                        richness_estimator = richness_estimator,
                        jack_alpha = 0.05, 
                        jack_max = 10,
                        coverage_estimator = coverage_estimator,
                        q = 0,
                        check_arguments = TRUE
                      )
                    )
                  }
                }
              ) 
              # Make a dataframe with the list to avoid nested lists
              the_df <- dplyr::bind_rows(the_list)
             }
          )
          # Make a dataframe with the list to avoid nested lists
          the_df <- dplyr::bind_rows(the_list)
          }
        )
        # Make a dataframe with the list to avoid nested lists
        the_df <- dplyr::bind_rows(the_list)
      }
    )
    # Coerce to a dataframe
    probabilities.dataframe <- dplyr::bind_rows(probabilities.list)
    
    # All probabilities must be below 1
    testthat::expect_lte(
      max(probabilities.dataframe$weight),
      1
    )
  }
)
