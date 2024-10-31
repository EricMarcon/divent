# Check all estimators

# Data: a vector of abundances
abd <- abd_species(paracou_6_abd[1, ])

testthat::test_that(
  "Naive estimator is correct", {
    testthat::skip_on_cran()
    # Naive
    testthat::expect_equal(
      div_richness(abd, estimator = "naive")$diversity,
      # Length of the vector
      sum(abd > 0)
    )
    # Jacknife
    testthat::expect_equal(
      div_richness(abd, estimator = "jackknife")$diversity,
      SPECIES::jackknife(as.matrix(abd_freq_count(abd)))$Nhat
    )
  }
)

# Jacknife
divent <- div_richness(abd, estimator = "jackknife")
species <- SPECIES::jackknife(as.matrix(abd_freq_count(abd)))
testthat::test_that(
  "Jacknife estimator is equal to that of SPECIES", {
    testthat::skip_on_cran()
    testthat::expect_equal(
      divent$diversity,
      species$Nhat
    )
  }
)

# Chao1 and iChao
ichao <- div_richness(abd, estimator = "iChao")
chao1 <- div_richness(abd, estimator = "Chao1")
testthat::test_that(
  "iChao estimator is greater than Chao1", 
  {
    testthat::skip_on_cran()
    testthat::expect_gt(
      ichao$diversity,
      chao1$diversity
    )
  }
)

# Rarefy
rarefy <- div_richness(paracou_6_abd, estimator = "rarefy")
testthat::test_that(
  "The rarefy estimator is reported correctly", 
  {
    testthat::skip_on_cran()
    testthat::expect_equal(
      rarefy$estimator,
      rep("rarefy", 4)
    )
  }
)

# Combine all parameters
abundances <- paracou_6_abd[1, ]

testthat::test_that(
  "No estimator fails", {
    testthat::skip_on_cran()
    # Estimate diversity systematically
    div_richness.list <- lapply(
      # All estimators
      eval(formals(divent:::div_richness.numeric)$estimator), 
      function(estimator) {
        the_list <-lapply(
          # All probability estimators
          eval(formals(divent:::div_richness.numeric)$probability_estimator), 
          function(probability_estimator) {
            the_list <- lapply(
              # All unveilings
              eval(formals(divent:::div_richness.numeric)$unveiling), 
              function(unveiling) {
                # Forbidden combination raises an error
                if ((estimator == "rarefy" & unveiling == "none")) {
                  NULL
                } else {
                  # print(paste(estimator, probability_estimator, unveiling))
                  suppressWarnings(
                    div_richness(
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
    div_richness.dataframe <- do.call(rbind, div_richness.list)
    
    # The min value must be the number of observed species
    testthat::expect_equal(
      min(div_richness.dataframe$diversity),
      sum(abundances[1, !colnames(abundances) %in% non_species_columns] > 0)
    )
  }
)


# Interpolation and extrapolation
sample_size <- abd_sum(abundances, as_numeric = TRUE)
levels <- c(sample_size / 2, round(sample_size * 1.5))

testthat::test_that(
  "No estimator fails during interpolation and extrapolation", {
    testthat::skip_on_cran()
    # Estimate diversity systematically
    div_richness.list <- lapply(
      # All estimators
      eval(formals(divent:::div_richness.numeric)$estimator), 
      function(estimator) {
        the_list <-lapply(
          # All probability estimators
          eval(formals(divent:::div_richness.numeric)$probability_estimator), 
          function(probability_estimator) {
            the_list <- lapply(
              # All unveilings
              eval(formals(divent:::div_richness.numeric)$unveiling), 
              function(unveiling) {
                the_list <-lapply(
                  # All levels
                  levels, 
                  function(level) {
                    # Forbidden combination raises an error
                    if ((estimator == "rarefy" & unveiling == "none")) {
                      NULL
                    } else {
                      # print(paste(estimator, probability_estimator, unveiling, level))
                      suppressWarnings(
                        div_richness(
                          abundances,
                          estimator = estimator,
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
    div_richness.dataframe <- do.call(rbind, div_richness.list)
    
    # The min value must be over 134 species
    testthat::expect_gt(
      min(div_richness.dataframe$diversity),
      134
    )
  }
)
