# Preparation ----
# Error message returned by the error_message() function
error_message_text <- "Check the function arguments."
# Function to simulate all argument checks.
# Must have the same arguments as check_divent_args, with a valid value
test_fun <- function(
  abd = NULL,
  abundances = paracou_6_abd,
  alpha = 0.05,
  as_numeric = TRUE,
  bootstrap = NULL,
  check_arguments = TRUE,
  correction = NULL,
  coverage_estimator = NULL,
  distances = diag(0, nrow = 2, ncol = 2),
  distribution = NULL,
  entropy_estimator = NULL,
  estimator = NULL,
  fisher_alpha = 40,
  gamma = TRUE,
  global = TRUE,
  jack_alpha = 0.05,
  jack_max = 10,
  k = 2,
  level = 1,
  n = 1,
  n_simulations = 2,
  normalize = TRUE,
  orders = 0,
  prob = rep(0.5, 2),
  prob_geom = 0.1,
  probability_estimator = NULL,
  q = 1,
  q_threshold = 10,
  r = 0:1,
  rate = 1,
  richness_estimator = NULL,
  sample_coverage = NULL,
  sd_lnorm = NULL,
  show_progress = TRUE,
  similarities = diag(nrow = 2, ncol = 2),
  size = NULL,
  species_distribution = paracou_6_abd,
  species_number = NULL,
  thomas_mu = NULL,
  thomas_scale = NULL,
  tree = NULL,
  use.names = NULL,
  unveiling = NULL,
  weights = NULL,
  w_max = 1,
  w_mean = 1,
  w_min = 1,
  weibull_scale = 30,
  weibull_shape = 2,
  X = paracou_6_wmppp,
  win = NULL
) {
  check_divent_args()
}

# check_divent_args() ----
testthat::test_that("Arguments are checked", {
  testthat::skip_on_cran()
  # abd
  testthat::expect_error(
    test_fun(abd = "a"),
    error_message_text
  )
  testthat::expect_error(
    test_fun(abd = -1),
    error_message_text
  )
  testthat::expect_error(
    test_fun(abd = 1 / 2),
    error_message_text
  )
  # abundances
  testthat::expect_error(
    test_fun(abundances = 0),
    error_message_text
  )
  # alpha
  testthat::expect_error(
    test_fun(alpha = "a"),
    error_message_text
  )
  testthat::expect_error(
    test_fun(alpha = -1),
    error_message_text
  )
  # as_numeric
  testthat::expect_error(
    test_fun(as_numeric = "a"),
    error_message_text
  )
  # check_arguments
  testthat::expect_error(
    test_fun(check_arguments = "a"),
    error_message_text
  )
  # distances
  testthat::expect_error(
    test_fun(distances = "a"),
    error_message_text
  )
  testthat::expect_error(
    test_fun(distances = 1),
    error_message_text
  )
  testthat::expect_error(
    test_fun(distances = matrix(0, ncol = 2, nrow = 3)),
    error_message_text
  )
  testthat::expect_error(
    test_fun(distances = matrix(c(0, 1, 1, 1), ncol = 2, nrow = 2)),
    error_message_text
  )
  testthat::expect_error(
    test_fun(distances = matrix(c(0, -1, -1, 0), ncol = 2, nrow = 2)),
    error_message_text
  )
  testthat::expect_error(
    test_fun(distances = matrix(1, ncol = 2, nrow = 2)),
    error_message_text
  )
  # fisher_alpha
  testthat::expect_error(
    test_fun(fisher_alpha = "a"),
    error_message_text
  )
  testthat::expect_error(
    test_fun(fisher_alpha = -1),
    error_message_text
  )
  # gamma
  testthat::expect_error(
    test_fun(gamma = "a"),
    error_message_text
  )
  # global
  testthat::expect_error(
    test_fun(global = "a"),
    error_message_text
  )
  # jack_alpha
  testthat::expect_error(
    test_fun(jack_alpha = "a"),
    error_message_text
  )
  testthat::expect_error(
    test_fun(jack_alpha = 2),
    error_message_text
  )
  # jack_max
  testthat::expect_error(
    test_fun(jack_max = "a"),
    error_message_text
  )
  testthat::expect_error(
    test_fun(jack_max = 20),
    error_message_text
  )
  # k
  testthat::expect_error(
    test_fun(k = "a"),
    error_message_text
  )
  testthat::expect_error(
    test_fun(k = -1),
    error_message_text
  )
  testthat::expect_error(
    test_fun(k = 1 / 2),
    error_message_text
  )
  # level
  testthat::expect_error(
    test_fun(level = "a"),
    error_message_text
  )
  testthat::expect_error(
    test_fun(level = -1),
    error_message_text
  )
  # n
  testthat::expect_error(
    test_fun(n = "a"),
    error_message_text
  )
  testthat::expect_error(
    test_fun(n = 0),
    error_message_text
  )
  # n_simulations
  testthat::expect_error(
    test_fun(n_simulations = "a"),
    error_message_text
  )
  # n_simulations
  testthat::expect_error(
    test_fun(n_simulations = 1),
    error_message_text
  )
  # normalize
  testthat::expect_error(
    test_fun(normalize = "a"),
    error_message_text
  )
  # orders
  testthat::expect_error(
    test_fun(orders = "a"),
    error_message_text
  )
  # prob
  testthat::expect_error(
    test_fun(prob = "a"),
    error_message_text
  )
  testthat::expect_error(
    test_fun(prob = -1),
    error_message_text
  )
  testthat::expect_error(
    test_fun(prob = 1:2),
    error_message_text
  )
  # prob_geom
  testthat::expect_error(
    test_fun(prob_geom = "a"),
    error_message_text
  )
  testthat::expect_error(
    test_fun(prob_geom = -1),
    error_message_text
  )
  testthat::expect_error(
    test_fun(prob_geom = 1:2),
    error_message_text
  )
  # q
  testthat::expect_error(
    test_fun(q = "a"),
    error_message_text
  )
  # q_threshold
  testthat::expect_error(
    test_fun(q_threshold = "a"),
    error_message_text
  )
  testthat::expect_error(
    test_fun(q_threshold = -1),
    error_message_text
  )
  testthat::expect_error(
    test_fun(q_threshold = 1 / 2),
    error_message_text
  )
  # r
  testthat::expect_error(
    test_fun(r = "a"),
    error_message_text
  )
  testthat::expect_error(
    test_fun(r = 1),
    error_message_text
  )
  testthat::expect_error(
    test_fun(r = 1:2),
    error_message_text
  )
  testthat::expect_error(
    test_fun(r = 0:-1),
    error_message_text
  )
  # rate
  testthat::expect_error(
    test_fun(rate = "a"),
    error_message_text
  )
  testthat::expect_error(
    test_fun(rate = -1),
    error_message_text
  )
  # sample_coverage
  testthat::expect_error(
    test_fun(sample_coverage = "a"),
    error_message_text
  )
  testthat::expect_error(
    test_fun(sample_coverage = 2),
    error_message_text
  )
  # sd_lnorm
  testthat::expect_error(
    test_fun(sd_lnorm = "a"),
    error_message_text
  )
  testthat::expect_error(
    test_fun(sd_lnorm = -1),
    error_message_text
  )
  # show_progress
  testthat::expect_error(
    test_fun(show_progress = "a"),
    error_message_text
  )
  # similarities
  testthat::expect_error(
    test_fun(similarities = "a"),
    error_message_text
  )
  testthat::expect_error(
    test_fun(similarities = 1),
    error_message_text
  )
  testthat::expect_error(
    test_fun(similarities = matrix(0, ncol = 2, nrow = 3)),
    error_message_text
  )
  testthat::expect_error(
    test_fun(similarities = matrix(c(0, 1, 1, 1), ncol = 2, nrow = 2)),
    error_message_text
  )
  testthat::expect_error(
    test_fun(similarities = matrix(c(0, -1, -1, 0), ncol = 2, nrow = 2)),
    error_message_text
  )
  testthat::expect_error(
    test_fun(similarities = matrix(0, ncol = 2, nrow = 2)),
    error_message_text
  )
  # size
  testthat::expect_error(
    test_fun(size = "a"),
    error_message_text
  )
  testthat::expect_error(
    test_fun(size = -1),
    error_message_text
  )
  # species_distribution
  testthat::expect_error(
    test_fun(species_distribution = "a"),
    error_message_text
  )
  # species_number
  testthat::expect_error(
    test_fun(species_number = "a"),
    error_message_text
  )
  testthat::expect_error(
    test_fun(species_number = -1),
    error_message_text
  )
  # thomas_mu
  testthat::expect_error(
    test_fun(thomas_mu = "a"),
    error_message_text
  )
  testthat::expect_error(
    test_fun(thomas_mu = -1),
    error_message_text
  )
  # thomas_scale
  testthat::expect_error(
    test_fun(thomas_scale = "a"),
    error_message_text
  )
  testthat::expect_error(
    test_fun(thomas_scale = -1),
    error_message_text
  )
  # tree
  testthat::expect_error(
    test_fun(tree = "a"),
    error_message_text
  )
  # weights
  testthat::expect_error(
    test_fun(weights = "a"),
    error_message_text
  )
  testthat::expect_error(
    test_fun(weights = -1),
    error_message_text
  )
  # w_max
  testthat::expect_error(
    test_fun(w_max = "a"),
    error_message_text
  )
  testthat::expect_error(
    test_fun(w_max = -1),
    error_message_text
  )
  # w_mean
  testthat::expect_error(
    test_fun(w_mean = "a"),
    error_message_text
  )
  testthat::expect_error(
    test_fun(w_mean = -1),
    error_message_text
  )
  # w_min
  testthat::expect_error(
    test_fun(w_min = "a"),
    error_message_text
  )
  testthat::expect_error(
    test_fun(w_min = -1),
    error_message_text
  )
  # weibull_scale
  testthat::expect_error(
    test_fun(weibull_scale = "a"),
    error_message_text
  )
  testthat::expect_error(
    test_fun(weibull_scale = -1),
    error_message_text
  )
  # weibull_shape
  testthat::expect_error(
    test_fun(weibull_shape = "a"),
    error_message_text
  )
  testthat::expect_error(
    test_fun(weibull_shape = -1),
    error_message_text
  )
  # X
  testthat::expect_error(
    test_fun(X = "a"),
    error_message_text
  )
  # win
  testthat::expect_error(
    test_fun(win = "a"),
    error_message_text
  )
})

# error_message() ----
testthat::test_that("error_message() works", {
  testthat::skip_on_cran()
  testthat::expect_error(
    error_message(
      message = "",
      argument = "",
      parent_function = "FUN"
    ),
    error_message_text
  )
})
