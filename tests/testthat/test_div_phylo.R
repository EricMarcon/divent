# Combine all parameters
abundances <- paracou_6_abd[1, ]
# Tree
tree <- paracou_6_taxo
# integer and non-integer q's
orders <- (0:6) / 2

testthat::test_that(
  "No estimator fails", {
    # Estimate diversity systematically
    div_phylo.list <- lapply(
      # All estimators
      eval(formals(divent:::div_phylo.numeric)$estimator), 
      function(estimator) {
        the_list <-lapply(
          # All probability estimators
          eval(formals(divent:::div_phylo.numeric)$probability_estimator), 
          function(probability_estimator) {
            the_list <-lapply(
              # All richness estimators
              eval(formals(divent:::div_phylo.numeric)$richness_estimator), 
              function(richness_estimator) {
                the_list <-lapply(
                  # All q's
                  orders, 
                  function(q) {
                    the_list <- lapply(
                      # All unveilings
                      eval(formals(divent:::div_phylo.numeric)$unveiling), 
                      function(unveiling) {
                        # print(paste(estimator, probability_estimator, unveiling, richness_estimator, q))
                        suppressWarnings(
                          div_phylo(
                            abundances,
                            tree = tree,
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
    div_phylo.dataframe <- do.call(rbind, div_phylo.list)
    
    # The min value must be UnveilJ / Chao2013 without unveiling
    testthat::expect_equal(
      min(div_phylo.dataframe$diversity, na.rm = TRUE),
      div_phylo(
        abundances, 
        tree = tree,
        q = max(orders), 
        probability_estimator = "Chao2013",
        unveiling = "none"
      )$diversity
    )
  }
)

# Interpolation and extrapolation
sample_size <- abd_sum(abundances, as_numeric = TRUE)
levels <- c(0.7, round(sample_size * 1.5))

testthat::test_that(
  "No estimator fails during interpolation and extrapolation", {
    # Estimate diversity systematically
    div_phylo.list <- lapply(
      # All probability estimators
      eval(formals(divent:::div_phylo.numeric)$probability_estimator), 
      function(probability_estimator) {
        the_list <-lapply(
          # All richness estimators
          eval(formals(divent:::div_phylo.numeric)$richness_estimator), 
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
                      eval(formals(divent:::div_phylo.numeric)$unveiling), 
                      function(unveiling) {
                        # print(paste(probability_estimator, unveiling, richness_estimator, q, level))
                        suppressWarnings(
                          div_phylo(
                            abundances,
                            tree = tree,
                            q = q,
                            # Estimator is not used
                            estimator = "naive",
                            level = level,
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
    div_phylo.dataframe <- do.call(rbind, div_phylo.list)
    
    # The min value must be UnveilJ / Chao2013 without unveiling
    testthat::expect_equal(
      min(div_phylo.dataframe$diversity, na.rm = TRUE),
      div_phylo(
        abundances, 
        tree = tree,
        q = max(orders), 
        level = min(levels)
      )$diversity
    )
  }
)
