# Combine all parameters
abundances <- paracou_6_abd[1, ]

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
            the_list <- lapply(
              # All unveilings
              eval(formals(divent:::div_hill.numeric)$unveiling), 
              function(unveiling) {
                print(paste(estimator, probability_estimator, unveiling))
                div_hill(
                  abundances,
                  q = 1.5,
                  estimator = estimator,
                  level = NULL,
                  probability_estimator = probability_estimator,
                  unveiling = unveiling,
                  as_numeric = FALSE,
                  check_arguments = TRUE
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
    div_hill.dataframe <- do.call(rbind, div_hill.list)
  }
)
