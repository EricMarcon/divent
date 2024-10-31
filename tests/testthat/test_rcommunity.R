# Combine all parameters
testthat::test_that(
  "No estimator fails", {
    testthat::skip_on_cran()
    # Simulate communities systematically
    rcommunity.list <- lapply(
      # All distributions
      eval(formals(divent:::rcommunity)$distribution), 
      function(distribution) {
        # print(distribution)
        rcommunity(
          1,
          size = 1000,
          species_number = 300,
          distribution = distribution,
          check_arguments = TRUE
        )
      }
    ) 

    # Coerce to a dataframe
    rcommunity.dataframe <- dplyr::bind_rows(rcommunity.list)
    # Replace NA's due to binding by zeros
    rcommunity.dataframe <- dplyr::mutate(
      rcommunity.dataframe,
      dplyr::across(
        .cols = dplyr::everything(),
        .fns = ~ ifelse(is.na(.x), 0, .x))
    )
    
    
    # The number of species must be less than 300
    testthat::expect_gte(
      300,
      ncol(abd_species(rcommunity.dataframe)),
    )
  }
)

# Unveil probabilities
abd <- abd_species(paracou_6_abd[1, ])

testthat::test_that(
  "No estimator fails", {
    testthat::skip_on_cran()
    # Simulate communities systematically
    rcommunity.list <- lapply(
      # All distributions
      eval(formals(divent:::rcommunity)$bootstrap), 
      function(bootstrap) {
        # print(bootstrap)
        rcommunity(
          1,
          abd = abd,
          bootstrap = bootstrap,
          check_arguments = TRUE
        )
      }
    ) 
    
    # Coerce to a dataframe
    rcommunity.dataframe <- dplyr::bind_rows(rcommunity.list)
    # Replace NA's due to binding by zeros
    rcommunity.dataframe <- dplyr::mutate(
      rcommunity.dataframe,
      dplyr::across(
        .cols = dplyr::everything(),
        .fns = ~ ifelse(is.na(.x), 0, .x))
    )
    
    
    # The number of individuals must equal the sample size
    testthat::expect_gte(
      unique(rcommunity.dataframe$weight),
      sum(abd)
    )
  }
)
