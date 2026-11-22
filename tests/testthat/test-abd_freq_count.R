# Abundance vector
the_abd <- as.numeric(paracou_6_abd[1, ])

test_that("Interpolation to sample size and coverage are equal", {
  testthat::skip_on_cran()
  # Coverage
  the_coverage <- 0.5
  # Size corresponding to coverage
  the_size <- coverage_to_size(
    the_abd,
    sample_coverage = the_coverage,
    as_numeric = TRUE
  )
  expect_identical(
    abd_freq_count(the_abd, level = the_coverage),
    abd_freq_count(the_abd, level = the_size)
  )
})

test_that("Extrapolation produces the correct number of singletons", {
  testthat::skip_on_cran()
  # Double the size of Paracou subplot 1
  the_size <- abd_sum(paracou_6_abd[1, ], as_numeric = TRUE)
  expect_equal(
    as.numeric(
      abd_freq_count(
        the_abd,
        level = the_size
        # Number of singletons
      )[1, "number_of_species"]),
    # Correct value
    84
  )
})
