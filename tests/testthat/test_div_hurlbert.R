# Combine all parameters
abundances <- paracou_6_abd[1, ]

testthat::test_that(
  "No estimator fails", {
    testthat::skip_on_cran()
    # Estimate diversity systematically
    div_hurlbert.list <- lapply(
      # All estimators
      eval(formals(divent:::div_hurlbert.numeric)$estimator), 
      FUN = function(estimator) {
        suppressWarnings(
          div_hurlbert(
            abundances,
            estimator = estimator,
            as_numeric = FALSE,
            check_arguments = TRUE
          )
        )
      }
    )
    # Coerce to a dataframe
    div_hurlbert.dataframe <- do.call(rbind, div_hurlbert.list)
    
    # The min value must be naive
    testthat::expect_equal(
      min(div_hurlbert.dataframe$diversity),
      div_hurlbert(
        abundances, 
        estimator = "naive"
      )$diversity
    )
  }
)
