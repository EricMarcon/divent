# Preparation ----
abd <- 1:2
abd_char <- c("a", "b", "b")

testthat::test_that("Objects are created correctly", {
  testthat::expect_equal(
    sum(as_probabilities(paracou_6_wmppp)[, -(1:2)]),
    1
  )
  testthat::expect_equal(
    sum(as_probabilities(abd_char)[, -(1:2)]),
    1
  )
  testthat::expect_equal(
    sum(as_probabilities(as.factor(abd_char))[, -(1:2)]),
    1
  )
  testthat::expect_equal(
    sum(as_abundances(paracou_6_wmppp)[, -(1:2)]),
    nrow(paracou_6_wmppp$marks)
  )
  testthat::expect_equal(
    sum(as_abundances(abd_char)[, -(1:2)]),
    length(abd_char)
  )
  testthat::expect_equal(
    sum(as_abundances(as.factor(abd_char))[, -(1:2)]),
    length(abd_char)
  )
  testthat::expect_equal(
    abd_sum(abundances(abd), as_numeric = TRUE),
    sum(abd)
  )
  testthat::expect_equal(
    sum(prob_species(as_probabilities(1))),
    1
  )
  testthat::expect_equal(
    nrow(prob_species(as_probabilities(diag(nrow = 2)))),
    2
  )
  testthat::expect_equal(
    sum(as.numeric(as_species_distribution(abd)) - abd),
    0
  )
  testthat::expect_equal(
    sum(as_named_vector.wmppp(paracou_6_wmppp)),
    abd_sum(as_abundances(paracou_6_wmppp), as_numeric = TRUE)
  )
  testthat::expect_error(
    as_probabilities(-1),
    "Species probabilities must be positive."
  )

  testthat::expect_error(
    abundances("a"),
    "`x` must be numeric"
  )
  testthat::expect_error(
    abundances(-1),
    "Species abundances must be positive."
  )
  testthat::expect_condition(
    as.numeric(as_species_distribution(as_probabilities(diag(nrow = 2)))),
    "The species_distribution object contains several rows."
  )
})
