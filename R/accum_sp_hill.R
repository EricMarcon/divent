#' Spatial Diversity Accumulation of a Community
#'
#' Spatial Diversity and Entropy Accumulation Curves represent the accumulation of
#' entropy and diversity with respect to the distance from individuals
#'
#' `accum_sp_hill()` or `accum_sp_tsallis()` estimate the diversity or entropy
#' accumulation curve of a distribution.
#'
#' @inheritParams check_divent_args
#'
#' @name accum_sp_hill
NULL


#' @rdname accum_sp_hill
#'
#' @export
#' @param orders A numeric vector: the diversity orders to address. Default is 0.
#' @param neighbors A vector of integers.
#' Entropy will be accumulated along this number of neighbors around each individual.
#' Default is 10% of the individuals.
#' @param r A vector of distances.
#' If `NULL` accumulation is along `n`, else neighbors are accumulated in circles of radius `r`.
#' @param correction The edge-effect correction to apply when estimating
#' the entropy of a neighborhood community that does not fit in the window.
#' Does not apply if neighborhoods are defined by the number of neighbors.
#' Default is "none".
#' "extrapolation" extrapolates the observed diversity up to the number of individuals
#' estimated in the full area of the neighborhood, which is slow.
#' @param individual If `TRUE`, individual neighborhood entropies are returned.
#'
#' @return An [accum_sp] object, that is also either an [accum_sp_diversity],
#' [accum_sp_entropy] or [accum_sp_mixing] object.
#'
#' @export
#'
#' @examples
#' # Generate a random community
#' X <- rspcommunity(1, size = 50, species_number = 3)
#' # Calculate the accumulation of richness
#' accum <- accum_sp_hill(X)
#' plot(accum, q = 0)
#' # along distance
#' accum_r <- accum_sp_hill(X, orders = 1, r = seq(0, .5, .05))
#' autoplot(accum_r, q = 1)
#'
accum_sp_tsallis <- function(
    X,
    orders = 0,
    neighbors = 1:ceiling(X$n / 2),
    r = NULL,
    correction = c("none", "extrapolation"),
    richness_estimator = c("rarefy", "jackknife", "iChao1", "Chao1", "naive"),
    individual = FALSE,
    show_progress = TRUE,
    check_arguments = TRUE) {

  # Check arguments
  richness_estimator <- match.arg(richness_estimator)
  correction <- match.arg(correction)
  if (any(check_arguments)) {
    check_divent_args()
  }

  if (is.null(r)) {
    # A number of neighbors ----
    # n nearest neighbors. Find them.
    neighbors.matrix <- spatstat.geom::nnwhich(X, k = neighbors)
    # Add the reference point to get a table: center points in rows, neighbors in columns,
    # including the point itself in the first column
    neighbors.matrix <- cbind(Reference = 1:X$n, neighbors.matrix)

    # Prepare a progress bar and the result arrays
    if (show_progress & interactive()) {
      cli::cli_progress_bar("Computing entropy", total = length(neighbors))
    }
    # 3D array: q, n, observed entropy in a single slice
    ent_q_nr_observed <- array(
      0,
      dim = c(length(orders), 1 + length(neighbors), 1),
      dimnames = list(q = orders, n = c(1, 1 + neighbors), Values = "Observed")
    )
    # individual values. 3D array: q, n, individuals
    if (individual) {
      ent_q_nr_individuals <- array(
        0,
        dim = c(length(orders), 1 + length(neighbors), X$n),
        dimnames = list(
          q = orders,
          n = c(1, 1 + neighbors),
          Point = row.names(X$marks)
        )
      )

    } else {
      ent_q_nr_individuals <- NA
    }

    # At each number of neighbors, calculate the entropy of
    # all points' neighborhood for each q
    for (k in seq_along(neighbors)) {
      # Neighbor communities: 1 community per column
      neighbor_communities <- apply(
        neighbors.matrix[, 1:(k + 1)],
        MARGIN = 1,
        FUN = function(neighbors) spatstat.geom::marks(X)$PointType[neighbors]
      )
      # Calculate entropy of each neighborhood and all q values
      ent_nbhood_q <- apply(
        neighbor_communities,
        MARGIN = 2,
        FUN = function(community) {
          sapply(
            orders,
            function(q) {
              ent_tsallis(
                as_abundances.character(community),
                q = q,
                estimator = "naive",
                as_numeric = TRUE,
                check_arguments = FALSE
              )
            }
          )
        }
      )
      # Keep individual neighborhood values
      if (individual) {
        ent_q_nr_individuals[, k + 1, ] <- ent_nbhood_q
      }

      # Mean entropy.
      # If ent_nbhood_q is a vector (i.e. a single value of q is provided),
      # transpose it to get a 1-row matrix.
      if (is.null(dim(ent_nbhood_q))) {
        ent_nbhood_q <- t(ent_nbhood_q)
      }
      ent_q_nr_observed[, k + 1, 1] <- apply(
        t(t(ent_nbhood_q)),
        MARGIN = 1,
        FUN = mean,
        na.rm = TRUE
      )
      if (show_progress & interactive()) cli::cli_progress_update()
    }
    if (show_progress & interactive()) cli::cli_progress_done()

    # Entropy of a single individual is 0.
    # This is the default value of the arrays so don't run.
    #  ent_q_nr_observed[, 1, 1] <- 0
    #  if (individual) ent_q_nr_individuals[, 1, ] <- 0

  } else {
    # A vector of distances ----
    # neighbors up to distance r. Distances are in r.
    # The first distance is 0 (verified by check_arguments)

    # The max value of the factors is needed
    species_number <- max(as.integer(spatstat.geom::marks(X)$PointType))
    # Run C++ routine to fill a 3D array.
    # Rows are points, columns are r, the 3rd dimension has a z-value per species.
    # Values are the number (weights) of neighbors of each point,
    # up to distance r, of species z.
    neighbors.array <- parallelCountNbd(
      r = r,
      NbSpecies = species_number,
      x = X$x,
      y = X$y,
      Type = spatstat.geom::marks(X)$PointType,
      Weight = spatstat.geom::marks(X)$PointWeight
    )
    # The array of neighbor communities is built from the vector returned.
    dim(neighbors.array) <- c(X$n, length(r), species_number)

    # Prepare a progress bar and the result arrays
    if (show_progress & interactive()) {
      cli::cli_progress_bar("Computing entropy", total = length(r))
    }
    # 3D array: q, r, observed entropy in a single slice
    ent_q_nr_observed <- array(
      0,
      dim = c(length(orders), length(r), 1),
      dimnames = list(q = orders, r = r, Values = "Observed")
    )
    # individual values
    if (individual) {
      ent_q_nr_individuals <- array(
        0,
        dim = c(length(orders), length(r), X$n),
        dimnames = list(
          q = orders,
          r = r,
          Point = row.names(spatstat.geom::marks(X))
        )
      )
    } else {
      ent_q_nr_individuals <- NA
    }

    # At each distance, calculate the entropy of
    # all points' neighborhood for each q
    for (d in 2:length(r)) {
      # Neighbor communities of each point at distance r: 1 species per column
      neighbor_communities <- vapply(
        1:species_number,
        FUN = function(i) rowSums(neighbors.array[, 1:d, i]),
        FUN.VALUE = numeric(X$n)
      )
      # Calculate entropy of each community and all q values
      if (correction == "none") {
        # No edge-effect correction
        ent_nbhood_q <-  apply(
          neighbor_communities,
          MARGIN = 1,
          FUN = function(community) {
            vapply(
              orders,
              FUN = function(q) {
                ent_tsallis(
                  as_abundances.numeric(community),
                  q = q,
                  estimator = "naive",
                  as_numeric = TRUE,
                  check_arguments = FALSE
                )
              },
              FUN.VALUE = 0
            )
          }
        )
      } else {
        if (correction == "extrapolation") {
          # Number of neighbors of each point
          neighbors.matrix <- rowSums(neighbor_communities)
          # Edge effects
          the_extrapolation <- integer(X$n)
          for (i in 1:X$n) {
            # Intersection between the point's neighborhood and the window
            the_intersection <- spatstat.geom::area(
              spatstat.geom::intersect.owin(
                X$window,
                spatstat.geom::disc(radius = r[d], centre = c(X$x[i], X$y[i]))
              )
            )
            # the_extrapolation ratio is that of the whole disc
            # to the part of the disc inside the window
            the_extrapolation[i] <- as.integer(
              neighbors.matrix[i] * pi * r[d]^2 / the_intersection
            )
          }
          # Prepare an array to store the results
          ent_nbhood_q <- array(
            0,
            dim = c(length(orders), nrow(neighbor_communities))
          )
          for (community in 1:nrow(neighbor_communities)) {
            for (order in seq_along(orders)) {
              # Suppress the warnings for Coverage=0 every time neighbors are singletons only.
              suppressWarnings(
                ent_nbhood_q[order, community] <- ent_tsallis(
                  neighbor_communities[community, ],
                  q = orders[order],
                  level = the_extrapolation[community],
                  richness_estimator = richness_estimator, # TODO: check estimator
                  as_numeric = TRUE,
                  check_arguments = FALSE
                )
              )
            }
          }
        } else {
          cli::cli_abort(
            "The edge-effect correction argument correction has not been recognized."
          )
        }
      }
      # Keep individual neighborhood values
      if (individual) {
        ent_q_nr_individuals[, d, ] <- ent_nbhood_q
      }
      # Mean entropy.
      # If ent_nbhood_q is a vector (i.e. a single value of q is provided),
      # transpose it to get a 1-row matrix.
      if (is.null(dim(ent_nbhood_q))) {
        ent_nbhood_q <- t(ent_nbhood_q)
      }
      ent_q_nr_observed[, d, 1] <- apply(
        t(t(ent_nbhood_q)),
        MARGIN = 1,
        FUN = mean,
        na.rm = TRUE
      )
      if (show_progress & interactive()) cli::cli_progress_update()
    }
    if (show_progress & interactive()) cli::cli_progress_done()
    # Entropy at r=0 is 0. This is the default value of the arrays so don't run.
    #  ent_q_nr_observed[, 1, 1] <- 0
    #  if (individual) ent_q_nr_individuals[, 1, ] <- 0
  }

  entAccum <- list(
    X = X,
    accumulation = ent_q_nr_observed,
    neighborhoods = ent_q_nr_individuals
  )
  class(entAccum) <- c("accum_sp_entropy", "accum_sp", class(entAccum))
  return(entAccum)
}


#' @rdname accum_sp_hill
#' @param h0 The null hypothesis to compare the distribution of `X` to.
#' If "none", the default value, no null hypothesis is tested.
#' "multinomial" means the community will be rarefied down to the number of `neighbors`.
#' "random location" means the points will we randomly permuted across their actual locations.
#' "binomial" means the points will we uniformly and independently drawn
#' in the window (a binomial point process is a Poisson point process conditionally to the number of points).
#'
#' @export
#'
accum_sp_hill <- function(
    X,
    orders = 0,
    neighbors = 1:ceiling(X$n / 2),
    r = NULL,
    correction = c("none", "extrapolation"),
    richness_estimator = c("rarefy", "jackknife", "iChao1", "Chao1", "naive"),
    h0 = c("none", "multinomial", "random location", "binomial"),
    alpha = 0.05,
    n_simulations = 100,
    individual = FALSE,
    show_progress = TRUE,
    check_arguments = TRUE) {

  # Check arguments
  richness_estimator <- match.arg(richness_estimator)
  correction <- match.arg(correction)
  h0 <- match.arg(h0)
  if (any(check_arguments)) {
    check_divent_args()
  }

  # Prepare an array to store data
  if (is.null(r)) {
    # Neighborhoods defined as the number of neighbors + a column for no neighbor.
    n_cols <- 1 + length(neighbors)
    the_seq <- c(1, neighbors + 1)
  } else {
    # Neighborhoods defined by distances. The first distance in r is 0.
    n_cols <- length(r)
    the_seq <- r
  }

  # Get the entropy
  the_diversity <- accum_sp_tsallis(
    X = X,
    orders = orders,
    neighbors = neighbors,
    r = r,
    correction = correction,
    richness_estimator = richness_estimator,
    individual = individual,
    show_progress = (show_progress & (h0 == "none" | h0 == "multinomial")),
    check_arguments = FALSE)

  if (h0 == "none") {
    is_h0_found <- TRUE
  } else {
    # H0 will have to be found
    is_h0_found <- FALSE
    # Rename accumulation
    names(the_diversity)[2] <- "Entropy"
    # Put the entropy into a 4-D array.
    # 4 z-values: observed, expected under H0, lower and upper bounds of H0.
    the_diversity$accumulation <- rep(the_diversity$Entropy, 4)
    dim(the_diversity$accumulation) <- c(length(orders), n_cols, 4)
    dimnames(the_diversity$accumulation) <- list(
      q = orders,
      n = the_seq,
      c("Observed", "Theoretical", "Lower bound", "Upper bound")
    )
    # if accumulation is along r, change the name
    if (!is.null(r)) {
      names(dimnames(the_diversity$accumulation))[2] <- "r"
    }
    the_diversity$Entropy <- NULL
  }

  # Calculate Hill Numbers, by row
  for (order in seq_along(orders)) {
    # Transform entropy to diversity, by row (where q does not change)
    the_diversity$accumulation[order, , 1] <- exp_q(
      the_diversity$accumulation[order, , 1],
      q = orders[order]
    )
    if (individual) {
      the_diversity$neighborhoods[order, , ] <- exp_q(
        the_diversity$neighborhoods[order, , ],
        q = orders[order]
      )
    }
  }

  # Null distribution
  if (h0 == "multinomial") {
    # Rarefy the community
    if (!is.null(r)) {
      cli::cli_abort(
        paste(
          "The 'multinomial' null hypothesis only applies to accumulation",
          "by number of neighbors."
        )
      )
    }
    is_h0_found <- TRUE
    # Prepare a progress bar
    if (show_progress & interactive()) {
      cli::cli_progress_bar("Running simulations", total = length(orders))
    }

    # Prepare the distribution of the abundances of species.
    abd <- as_abundances.wmppp(X)
    for (order in seq_along(orders)) {
      # Rarefy the community to the sizes of neighborhoods
      h0_values <- accum_hill(
        abd,
        q = as.numeric(orders[order]),
        levels = the_seq,              # TODO: missing arguments for accum_hill.numeric
        n_simulations = n_simulations,
        alpha = alpha,
        show_progress = FALSE,
        check_arguments = FALSE
      )
      # Extract the results from the object returned
      the_diversity$accumulation[order, , 2] <- h0_values$diversity
      the_diversity$accumulation[order, , 3] <- h0_values$inf
      the_diversity$accumulation[order, , 4] <- h0_values$sup
      if (show_progress & interactive()) cli::cli_progress_update()
    }
    if (show_progress & interactive()) cli::cli_progress_done()
  }
  if (h0 == "random location" | h0 == "binomial") {
    is_h0_found <- TRUE
    # Prepare a progress bar
    if (show_progress & interactive()) {
      cli::cli_progress_bar("Running simulations", total = n_simulations)
    }
    # Prepare a 3-D array to store results. Rows are q, columns are r or n,
    # z-values are for each simulation.
    h0_diversity <- array(
      0,
      dim = c(length(orders), n_cols, n_simulations)
    )
    # Simulate communities according to H0
    for (i in (1:n_simulations)) {
      # Random community
      if (h0 == "random location")
        h0_X <- dbmss::rRandomLocation(X, CheckArguments = FALSE)
      if (h0 == "binomial")
        h0_X <- dbmss::rRandomPositionK(X, CheckArguments = FALSE)
      # Calculate its accumulated diversity
      h0_diversity[, , i] <- accum_sp_hill(
        h0_X,
        orders = orders,
        neighbors = neighbors,
        r = r,
        correction = correction,
        richness_estimator = richness_estimator,
        h0 = "none",
        n_simulations = 0,                 # TODO: check argument list
        individual = FALSE,
        show_progress = FALSE,
        check_arguments = FALSE
      )$accumulation[, , 1]
      if (show_progress & interactive()) cli::cli_progress_update()
    }
    if (show_progress & interactive()) cli::cli_progress_done()
    # Calculate quantiles
    for (q in seq_along(orders)) {
      for (r in seq_along(r)) {
        the_diversity$accumulation[q, r, 3:4] <- stats::quantile(
          h0_diversity[q, r, ], c(alpha, 1 - alpha), na.rm = TRUE
        )
        the_diversity$accumulation[q, r, 2] <- mean(
          h0_diversity[q, r, ], na.rm = TRUE
        )
      }
    }
  }

  if (!is_h0_found) {
    cli::cli_abort(
      "The value of 'h0' does not correspond to a valid null hypothesis."
    )
  }

  class(the_diversity) <- c("accum_sp_diversity", "accum_sp", class(the_diversity))
  return(the_diversity)
}


#' @rdname accum_sp_hill
#'
#' @export
#'
accum_mixing <- function(
    X,
    orders = 0,
    neighbors = 1:ceiling(X$n / 2),
    r = NULL,
    correction = c("none", "extrapolation"),
    richness_estimator = c("rarefy", "jackknife", "iChao1", "Chao1", "naive"),
    h0 = c("none", "multinomial", "random location", "binomial"),
    alpha = 0.05,
    n_simulations = 100,
    individual = FALSE,
    show_progress = TRUE,
    check_arguments = TRUE) {

  # Check arguments
  richness_estimator <- match.arg(richness_estimator)
  correction <- match.arg(correction)
  h0 <- match.arg(h0)
  if (any(check_arguments)) {
    check_divent_args()
  }

  # Get the diversity accumulation
  the_mixing <- accum_sp_hill(
    X,
    orders = orders,
    neighbors = neighbors,
    r = r,
    correction = correction,
    richness_estimator = richness_estimator,
    h0 = h0,
    n_simulations = n_simulations,                 # TODO: check argument list
    individual = show_progress,
    show_progress = show_progress,
    check_arguments = FALSE
  )

  # Normalize it
  the_mixing$accumulation[, , 1] <- the_mixing$accumulation[, , 1] /
    the_mixing$accumulation[, , 2]
  the_mixing$accumulation[, , 3] <- the_mixing$accumulation[, , 3] /
    the_mixing$accumulation[, , 2]
  the_mixing$accumulation[, , 4] <- the_mixing$accumulation[, , 4] /
    the_mixing$accumulation[, , 2]
  # Normalize individual values
  if (individual) {
    for (i in seq_len(X$n)) {
      the_mixing$neighborhoods[, , i] <- the_mixing$neighborhoods[, , i] /
        the_mixing$accumulation[, , 2]
    }
  }
  the_mixing$accumulation[, , 2] <- 1

  class(the_mixing) <- c("accum_sp_mixing", "accum_sp", class(the_mixing))
  return(the_mixing)
}
