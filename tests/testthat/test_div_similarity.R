# Combine all parameters
abundances <- paracou_6_abd[1, ]
# Similarities
Z <- fun_similarity(paracou_6_fundist)
# integer and non-integer q's
orders <- (0:6) / 2

testthat::test_that(
  "No estimator fails", {
    testthat::skip_on_cran()
    # Estimate diversity systematically
    div_similarity.list <- lapply(
      # All estimators
      eval(formals(divent:::div_similarity.numeric)$estimator), 
      function(estimator) {
        the_list <-lapply(
          # All probability estimators
          eval(formals(divent:::div_similarity.numeric)$probability_estimator), 
          function(probability_estimator) {
            the_list <-lapply(
              # All q's
              orders, 
              function(q) {
                the_list <- lapply(
                  # All unveilings
                  eval(formals(divent:::div_similarity.numeric)$unveiling), 
                  function(unveiling) {
                    # print(paste(estimator, probability_estimator, unveiling, richness_estimator, q))
                    suppressWarnings(
                      div_similarity(
                        abundances,
                        similarities = Z,
                        q = q,
                        estimator = estimator,
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
        # Make a dataframe with the list to avoid nested lists
        the_df <- do.call(rbind, the_list)
      }
    )
    # Coerce to a dataframe
    div_similarity.dataframe <- do.call(rbind, div_similarity.list)
    
    # The min value must be UnveilJ / Chao2013 without unveiling
    testthat::expect_equal(
      min(div_similarity.dataframe$diversity, na.rm = TRUE),
      div_similarity(
        abundances, 
        similarities = Z,
        q = 0.5, 
        probability_estimator = "Chao2013",
        unveiling = "none"
      )$diversity
    )
  }
)
