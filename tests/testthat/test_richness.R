# Check all estimators

# Data: a vector of abundances
abd <- abd_species(paracou_6_abd[1, ])

testthat::test_that(
  "Naive estimator is correct", 
  {
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
  "Jacknife estimator is equal to that of SPECIES", 
  {
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

# rarefy
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
