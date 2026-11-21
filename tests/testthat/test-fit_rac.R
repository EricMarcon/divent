# Draw random distribution ----
## Lognormal ----
sd_lnorm <- 2
lnorm_dist <- rcommunity(
  n = 1,
  size = 1E5,
  species_number = 200,
  distribution = "lnorm",
  sd_lnorm = sd_lnorm
)
# Numeric output
lnorm_dist_num <- as.numeric(lnorm_dist)
## Log-series ----
fisher_alpha <- 15
lseries_dist <- rcommunity(
  n = 1,
  size = 1E6,
  distribution = "lseries",
  fisher_alpha = fisher_alpha
)
# Numeric output
lseries_dist_num <- as.numeric(lseries_dist)
## Geometric ----
prob_geom  <- 0.05
geom_dist <- rcommunity(
  n = 1,
  size = 1E5,
  distribution = "geom",
  prob_geom  = prob_geom,
)
# Numeric output
geom_dist_num <- as.numeric(geom_dist)
## Broken stick ----
bstick_dist <- rcommunity(
  n = 1,
  size = 1E5,
  species_number = 200,
  distribution = "bstick"
)
# Numeric output
bstick_dist_num <- as.numeric(bstick_dist)

# Infer the parameters ----
lnorm_param <- fit_rac(lnorm_dist, distribution = "lnorm")
lnorm_num_param <- fit_rac(lnorm_dist, distribution = "lnorm")
lseries_param <- fit_rac(lseries_dist, distribution = "lseries")
lseries_num_param <- fit_rac(lseries_dist, distribution = "lseries")
geom_param <- fit_rac(geom_dist, distribution = "geom")
geom_num_param <- fit_rac(geom_dist, distribution = "geom")
bstick_param <- fit_rac(bstick_dist, distribution = "bstick")
bstick_num_param <- fit_rac(bstick_dist, distribution = "bstick")

# Tests ----
testthat::test_that("Distribution parameters are inferred correctly", {
  testthat::skip_on_cran()
  ## Lognormal ----
  testthat::expect_equal(
    lnorm_param$sigma,
    lnorm_num_param$sigma
  )
  testthat::expect_lt(
    abs(lnorm_param$sigma - sd_lnorm),
    expected = sd_lnorm / 10
  )
  ## Log-series ----
  testthat::expect_equal(
    lseries_param$alpha,
    lseries_num_param$alpha
  )
  testthat::expect_lt(
    abs(lseries_param$alpha - fisher_alpha),
    expected = fisher_alpha / 10
  )
  ## Geometric ----
  testthat::expect_equal(
    geom_param$prob,
    geom_num_param$prob
  )
  testthat::expect_lt(
    abs(geom_param$prob - prob_geom),
    expected = prob_geom / 10
  )
  ## Broken stick ----
  testthat::expect_equal(
    bstick_param$max,
    bstick_num_param$max
  )
})

testthat::test_that("Negative probabilities raise an error", {
  testthat::skip_on_cran()
  testthat::expect_error(
    # Abundance distribution with negative value (not checked to avoid error)
    fit_rac(as_abundances(c(1:3, -1), check_arguments = "FALSE")),
    "Species probabilities or abundances must be positive."
  )
  testthat::expect_error(
    # Numeric distribution
    fit_rac(c(1:3, -1)),
    "Species probabilities or abundances must be positive."
  )
})
