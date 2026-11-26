# Abundance vector
abd <- 2:3

# Run simulations for coverage
profile_hill(abd, n_simulations = 2)

testthat::test_that("Arguments are checked", {
  testthat::skip_on_cran()
  testthat::expect_error(
    profile_hill(-abd),
    "Species probabilities or abundances must be positive."
  )
  testthat::expect_error(
    profile_hill(abd, as_numeric = TRUE, n_simulations = 10),
    "No simulations are allowed if a numeric vector is expected."
  )
})
