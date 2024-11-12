#' Spatial Diversity Accumulation of a Community
#' 
#' Spatial Diversity and Entropy Accumulation Curves represent the accumulation of 
#' entropy and diversity with respect to the distance from individuals
#' 
#' `accum_sp_hill()` or `accum_sp_tsallis()` estimate the diversity or entropy 
#' accumulation curve of a distribution.
#' 
#' @inheritParams check_divent_args
#' @references
#' \insertAllCited{}
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
#' @return An "Accumulation" object that is a list. 
#' 
#' - Its first item, named "X", is `X`. 
#' - Its second item, named "Accumulation", is a 3-D array containing average entropy.
#' The third dimension of the array is only of length 1: it contains observed entropy.
#' The first two dimensions are respectively for $q$ values and the number of points 
#' of the neighborhood, starting from 1 (the point itself, with no neighbor), or the distances starting from 0.
#' - Its third item, named "Neighborhoods" has the same structure as the second one 
#' but its third dimension contains the local values accumulated in the neighborhood of each point. 
#'
#' @export
#'
#' @examples
#' # Generate a random community
#' X <- rspcommunity(1, size = 50, species_number = 3)
#' # Calculate the accumulation of richness 
#' accum <- accum_sp_tsallis(X)
#' # plot(accum, q=0)
#' # along distance
#' accumR <- accum_sp_tsallis(X, orders = 1, r = seq(0, .5, .05))
#' # plot(accumR, q=1)
#' 
accum_sp_tsallis <- function(
    X, 
    orders = 0, 
    estimator = c("UnveilJ", "ChaoJost", "ChaoShen", "GenCov", "Grassberger", 
                  "Marcon", "UnveilC", "UnveiliC", "ZhangGrabchak", "naive",
                  "Bonachela", "Holste"),
    neighbors = 1:ceiling(X$n / 2), 
    r = NULL, 
    correction = c("none", "extrapolation"),
    richness_estimator = c("rarefy", "jackknife", "iChao1", "Chao1", "naive"),
    individual = FALSE, 
    show_progress = TRUE, 
    check_arguments = TRUE) {
  
  if (any(check_arguments)) {
    check_divent_args()
  }
  estimator <- match.arg(estimator) 
  richness_estimator <- match.arg(richness_estimator) 
  correction <- match.arg(correction) 
  
  if (is.null(r)) {
    # A number of neighbors ----
    # n nearest neighbors. Find them.
    neighbors.matrix <- spatstat.geom::nnwhich(X, k = neighbors)
    # Add the reference point to get a table: center points in line, neighbors in columns, 
    # including the point itself in the first column
    neighbors.matrix <- cbind(Reference = 1:X$n, neighbors.matrix)
    
    # Prepare a progress bar and the result arrays
    if (show_progress & interactive()) {
      pgb <- utils::txtProgressBar(min = 0, max = length(neighbors))
    }
    # 3D array: q, n, observed entropy in a single slice
    ent_q_n_observed <- array(
      0, 
      dim = c(length(orders), 1 + length(neighbors), 1),
      dimnames = list(q = orders, n = c(1, 1 + neighbors), Values = "Observed")
    )
    # Individual values. 3D array: q, n, individuals
    if (individual) {
      ent_q_n_individuals <- array(
        0, 
        dim = c(length(orders), 1 + length(neighbors), X$n),
        dimnames = list(
          q = orders, 
          n = c(1, 1 + neighbors), 
          Point = row.names(X$marks)
        )
      )
      
    } else {
      ent_q_n_individuals <- NA
    }

    # At each number of neighbors, calculate the entropy of 
    # all points' neighborhood for each q
    for (k in 1:length(neighbors)) {
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
                as_abundances(community), 
                q = q, 
                estimator = "naive", # TODO: check the estimator
                as_numeric = TRUE,
                check_arguments = FALSE
              )
            }
          )
        } 
      )
      # Keep individual neighborhood values
      if (individual) {
        ent_q_n_individuals[, k + 1, ] <- ent_nbhood_q
      }

      # Mean entropy. 
      # If ent_nbhood_q is a vector (i.e. a single value of q is provided), 
      # transpose it to get a 1-row matrix.
      if (is.null(dim(ent_nbhood_q))) {
        ent_nbhood_q <- t(ent_nbhood_q)
      }
      ent_q_n_observed[, k + 1, 1] <- apply(
        t(t(ent_nbhood_q)), 
        MARGIN = 1, 
        FUN = mean, 
        na.rm = TRUE
      )
      if (show_progress & interactive()) {
        utils::setTxtProgressBar(pgb, k)
      }
    }
    if (show_progress & interactive()) close(pgb)
   
    # Entropy of a single individual is 0. 
    # This is the default value of the arrays so don't run.
    #  ent_q_n_observed[, 1, 1] <- 0
    #  if (individual) ent_q_n_individuals[, 1, ] <- 0
    
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
      pgb <- utils::txtProgressBar(min = 1, max = length(r))
    }
    # 3D array: q, r, observed entropy in a single slice
    ent_q_r_observed <- array(
      0, 
      dim = c(length(orders), length(r), 1),
      dimnames = list(q = orders, r = r, Values = "Observed")
    )
    # Individual values
    if (individual) {
      ent_q_r_individuals <- array(
        0, 
        dim = c(length(orders), length(r), X$n), 
        dimnames = list(
          q = orders, 
          r = r, 
          Point = row.names(spatstat.geom::marks(X))
        )
      )
    } else {
      ent_q_n_individuals <- NA
    }

    # At each distance, calculate the entropy of 
    # all points' neighborhood for each q
    for (d in 2:length(r)) {
      # Neighbor communities of each point at distance r: 1 community per column
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
                  as_abundances(community), 
                  q = q, 
                  estimator = "naive", # TODO: check the estimator
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
            for (order in 1:length(orders)) {
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
          stop("The edge-effect correction argument correction has not been recognized.")
        }
      }
      # Keep individual neighborhood values
      if (individual) {
        ent_q_n_individuals[, r, ] <- ent_nbhood_q
      }
      # Mean entropy. 
      # If ent_nbhood_q is a vector (i.e. a single value of q is provided), 
      # transpose it to get a 1-row matrix.
      if (is.null(dim(ent_nbhood_q))) {
        ent_nbhood_q <- t(ent_nbhood_q)
      }
      ent_q_r_observed[, r, 1] <- apply(
        t(t(ent_nbhood_q)), 
        MARGIN = 1, 
        FUN = mean, 
        na.rm = TRUE
      )
      if (show_progress & interactive()) utils::setTxtProgressBar(pgb, d)
    }
    if (show_progress & interactive()) close(pgb)
    # Entropy at r=0 is 0. This is the default value of the arrays so don't run.
    #  ent_q_r_observed[, 1, 1] <- 0
    #  if (individual) ent_q_n_individuals[, 1, ] <- 0
  }
  
  entAccum <- list(
    X = X, 
    Accumulation = ent_q_r_observed, 
    Neighborhoods = ent_q_n_individuals
  )
  class(entAccum) <- c("EntAccum", "Accumulation")
  return(entAccum)
}
