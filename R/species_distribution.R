#' Species Distributions
#' 
#' A Species Distribution is a [tibble::tibble] containing species abundances or probabilities.
#' 
#' `species_distribution` objects include `abundances` and `probabilities` objects.
#' 
#' `as_species_distribution()`, `as_abundances()`  and `as_probabilities` format 
#' the numeric, matrix or dataframe `x` so that appropriate 
#' versions of community functions (generic methods such as [plot] or 
#' [div_richness]) are applied. 
#' Abundance values are rounded (by default) to the nearest integer.
#' 
#' `as_probabilities()` normalizes the vector `x` so that it sums to 1. It gives
#' the same output as `probabilities()` with `estimator = "naive"`.
#' 
#' TODO: These functions can be applied to data frames to calculate the joint diversity \insertCite{Gregorius2010}{divent}.
#'  
#' `species_distribution` objects objects can be plotted by `plot()` and `autoplot()`.
#'
#' @param x An object.
#' @param ... Additional arguments to be passed to [plot]. Unused elsewhere.
#' @param check_arguments If `TRUE`, the function arguments are verified.
#' Should be set to `FALSE` to save time when the arguments have been checked elsewhere.
#' 
#' @examples
#' # Paracou data
#' is_species_distribution(paracou_6_abd)
#' # Whittaker plot fitted by a log-normal distribution
#' autoplot(paracou_6_abd[1,], distribution = "lognormal")
#' @references
#' \insertAllCited{}
#' 
#' 
#' @name species_distribution
NULL


#  Species Distribution ----

#' @rdname species_distribution
#'
#' @param names The names of the species distributions.
#' @param weights The weights of the sites of the species distributions.
#' 
#' @export
species_distribution <- function(
    x, 
    names = NULL, 
    weights = NULL, 
    check_arguments = TRUE) {

  # Check the data ----
  if (check_arguments) check_divent_args()
  if (!is.numeric(x)) stop("'x' must be numeric")
  if (length(dim(x)) > 2) stop("'x' may be a vector or a matrix")
  if (any(x < 0)) stop("Species distribution abundances or probabilities must be positive.")
  
  # Build a tibble from the data ----
  if (is.vector(x)) {
    ## Single distribution ----
    if (is.null(names(x))) {
      ### Columns: add default species names such as sp_1 ----
      names(x) <- paste(
        "sp", 
        formatC(seq_along(x), width = ceiling(log10(length(x))), flag = "0"),
        sep = "_"
      )
    }
    if (length(names) != 1) {
      ### Rows: Add a site name ----
      names <- paste("site", round(stats::runif(1)*.Machine$integer.max), sep="_")
    }
    # Build a tibble
    distribution <- tibble::as_tibble_row(c(site = names, x))
    
  } else {
    ## Several distributions ----
    if (is.null(colnames(x))) {
      ### Columns: add default species names such as sp_1 ----
      colnames(x) <- paste(
        "sp", 
        formatC(seq_len(ncol(x)), width = ceiling(log10(ncol(x))), flag = "0"),
        sep = "_"
      )
    }
    # Build a tibble
    distribution <- tibble::as_tibble(x, rownames = "site")
    ### Rows: site names = names or matrix row names or default ----
    if (!is.null(names)) {
      # site = names if the size matches
      if (length(names) == nrow(x)) {
        distribution$site <- names
      } else {
        stop("The length of 'names' must match the number of lines of the data matrix.")
      }
    } else {
      # names is null...
      if (is.null(row.names(x))) {
        # ...and no row names: set default names such as site_1
        distribution$site <- paste(
          "site", 
          formatC(seq_len(nrow(x)), width = ceiling(log10(nrow(x))), flag = "0"),
          sep = "_"
        )
      }
    }
    ### Rows: site weights ----
    if (!is.null(weights)) {
      # site = weights if the size matches
      if (length(weights) == nrow(x)) {
        distribution <- tibble::add_column(
          distribution, 
          weight = rowSums(x),
          .after = "site"
        )
      } else {
        stop("The length of 'weights' must match the number of lines of the data matrix.")
      }
    } else {
      # Weights are the number of individuals
      distribution <- tibble::add_column(
        distribution, 
        weight = rowSums(x),
        .after = "site"
      )
    }
  }
  
  # Set the class and return ----
  class(distribution) <- c("species_distribution", class(distribution))
  return(distribution)
}


#' @rdname species_distribution
#'
#' @export
as_species_distribution <- function(x, ...) {
  UseMethod("as_species_distribution")
}


#' @rdname species_distribution
#'
#' @export
as_species_distribution.numeric <- function(
    x,
    ...,
    check_arguments = TRUE) {
  
  return(
    species_distribution(
      x, 
      ..., 
      check_arguments = check_arguments
    )
  )
}


#' @rdname species_distribution
#'
#' @export
as_species_distribution.matrix <- function(
    x, 
    names = NULL, 
    weights = NULL, 
    ...,
    check_arguments = TRUE) {
  
  return(
    species_distribution(
      x, 
      names = names, 
      weights = weights, 
      ..., 
      check_arguments = check_arguments
    )
  )
}


#' @rdname species_distribution
#'
#' @export
as_species_distribution.data.frame <- function(
    x, 
    ...,
    check_arguments = TRUE) {
  
  if (check_arguments) check_divent_args()
  # Check the data
  if (any(x < 0)) stop("All numeric values of the dataframe must be positive.")
  
  # Build a tibble
  distribution <- tibble::as_tibble(x)
  
  # The first column should be "site"
  if (!"site" %in% colnames(distribution)) {
    distribution <- tibble::add_column(
      distribution, 
      site = paste(
        "site", 
        formatC(
          seq_len(nrow(distribution)), 
          width = ceiling(log10(nrow(distribution))), 
          flag = "0"
        ),
        sep = "_"
      ),
      .before = 1
    )
  }
  
  # The second column should be "weight"
  if (!"weight" %in% colnames(distribution)) {
    distribution <- tibble::add_column(
      distribution, 
      weight = rowSums(distribution[colnames(distribution) != "site"]),
      .after = "site"
    )
  }

  # Set the class and return
  class(distribution) <- c("species_distribution", class(distribution))
  return(distribution)
}


#' @rdname species_distribution
#'
#' @export
is_species_distribution <- function(x) {
  inherits(x, "species_distribution")
}


#' @rdname species_distribution
#'
#' @export
plot.species_distribution <- function(
    x, 
    ..., 
    distribution = c("lnorm", "lseries", "geom", "bstick"),
    log = "y", 
    main = NULL, 
    xlab = "Rank", 
    ylab = NULL) {
  
  # Prepare ylab
  if (is.null(ylab)) {
    if (is_probabilities(x)) {
      ylab <- "Probability"
    } else {
      ylab <- "Abundance" 
    }
  }

  # Find the max number of species in all communities
  abundances <- x[, !(colnames(x) %in% c("site", "weight"))] 
  s_obs_max <- max(rowSums(abundances > 0))
  
  # Prepare the plot: X and Y ranges
  graphics::plot(
    x = 1:s_obs_max,
    y = seq(from = 1, to = max(abundances), length.out = s_obs_max),
    type = "n",
    log = log, 
    main = main, 
    xlab = xlab, 
    ylab = ylab, 
    axes = FALSE, 
    ...
  )
  # X axis ticks must start from 1
  graphics::axis(1, graphics::axTicks(1) + 1)
  graphics::axis(2)
  graphics::box()
  
  # Color palette
  cols <- RColorBrewer::brewer.pal(nrow(x), "Set2")
  
  # Loop in communities to build the plot
  for (community in seq_len(nrow(x))) {
    # Extract the abundances of the community
    abd <- x[community, !(colnames(x) %in% c("site", "weight"))]
    # Eliminate zeros and sort
    abd <- sort(abd[abd > 0], decreasing = TRUE)
    sample_size <- sum(abd)
    s_obs <- length(abd)
    
    # Draw the species abundances
    graphics::points(
      x = seq_len(s_obs), 
      y = abd, 
      col = cols[community]
    )
    
    # Draw the fitted models
    if (!is.null(distribution)) {
      rac_fitted <- fit_rac(
        abd, 
        distribution = distribution, 
        check_arguments = FALSE
      )
      graphics::lines(
        x = rac_fitted$rac$rank, 
        y = rac_fitted$rac$abundance, 
        col = cols[community]
      )
    }
  }
  
  # Legend if several communities
  if (nrow(x) > 1) {
    legend(
      "topright",
      inset = .02,
      legend = x$site,
      col = cols,
      lty = 1,
      pch = 1
    )
  }
}


autoplot.species_distribution <- function(
    object, 
    ..., 
    distribution = c("lnorm", "lseries", "geom", "bstick"),
    ylog = TRUE, 
    main = NULL, 
    xlab = "Rank", 
    ylab = NULL, 
    pch = ggplot2::GeomPoint$default_aes$shape,
    cex = ggplot2::GeomPoint$default_aes$size) {

  # Prepare ylab
  if (is.null(ylab)) {
    if (is_probabilities(object)) {
      ylab <- "Probability"
    } else {
      ylab <- "Abundance" 
    }
  }
  
  # Find the max number of species in all communities
  s_obs_max <- max(
    rowSums(
      object[, !(colnames(object) %in% c("site", "weight"))] > 0
    )
  )
  
  # Prepare the plot
  the_plot <- ggplot2::ggplot() +
    ggplot2::scale_x_continuous(
      # X-axis starts at 0.01 to avoid the 0 X-label.
      limits = c(0.01, s_obs_max), 
      expand = c(0, 0)
    ) +
    ggplot2::labs(title = main, x = xlab, y = ylab)
  
  # Prepare the data to plot
  the_data <- data.frame(site = NULL, rank = NULL, abd = NULL)
  if (!is.null(distribution)) {
    the_model <- data.frame(site = NULL, rank = NULL, abd = NULL)
  }
  
  # Loop in communities to build the plot
  for (community in seq_len(nrow(object))) {
    # Extract the abundances of the community
    abd <- object[community, !(colnames(object) %in% c("site", "weight"))]
    # Eliminate zeros and sort
    abd <- sort(abd[abd > 0], decreasing = TRUE)
    sample_size <- sum(abd)
    s_obs <- length(abd)
    
    # Transform data into dataframe
    the_data <- rbind(
      the_data, 
      data.frame(
        site = object[community, "site"], 
        rank = seq_len(s_obs), 
        abundance = abd
      )
    )
    
    # Fitted model
    if (!is.null(distribution)) {
      rac_fitted <- fit_rac(
        abd, 
        distribution = distribution, 
        check_arguments = FALSE
      )
      the_model <- rbind(
        the_model,
        tibble::tibble(
          object[community, "site"],
          rac_fitted$rac
        )
      )
    }
  }
  
  # Plot
  the_plot <- the_plot +
    ggplot2::geom_point(
      data = the_data, 
      mapping = ggplot2::aes(x = .data$rank, y = .data$abundance, col = .data$site),
      shape = pch,
      size = cex
    ) 

  # Fitted models
  if (!is.null(distribution)) {
    the_plot <- the_plot + 
      ggplot2::geom_line(
        data = the_model,
        mapping = ggplot2::aes(x = .data$rank, y = .data$abundance, col = .data$site)
      )
  }
  
  # Log Y-axis
  if (ylog) the_plot <- the_plot + ggplot2::scale_y_log10() 
  
  # No legend if single community
  if (nrow(object) == 1) {
    the_plot <- the_plot + ggplot2::theme(legend.position = "none")
  }
  
  return(the_plot)
}


#  Probabilities ----


#' @rdname species_distribution
#'
#' @export
as_probabilities <- function(x, ...) {
  UseMethod("as_probabilities")
}


#' @rdname species_distribution
#'
#' @export
as_probabilities.numeric <- function(
    x,
    ..., 
    check_arguments = TRUE) {

  if (any(x < 0)) stop("Species probabilities must be positive.")

  # Normalize to 1
  prob <- x / sum(x)
  probabilities <- as_species_distribution(
    prob, 
    ...,
    check_arguments = check_arguments
  )

  class(probabilities) <- c("probabilities", class(probabilities))
  return(probabilities)
}


#' @rdname species_distribution
#' 
#' @export
as_probabilities.matrix <- function(
    x,
    names = NULL, 
    weights = NULL, 
    ...,
    check_arguments = TRUE) {
  
  probabilities <- as_species_distribution.matrix(
    x, 
    names = names, 
    weights = weights,
    ..., 
    check_arguments = check_arguments
  )
  
  class(probabilities) <- c("probabilities", class(probabilities))
  return(probabilities)
}


#' @rdname species_distribution
#'
#' @export
as_probabilities.data.frame <- function(
    x, 
    ...,
    check_arguments = TRUE) {
  
  probabilities <- as_species_distribution.data.frame(
    x, 
    ..., 
    check_arguments = check_arguments
  )
  
  class(probabilities) <- c("probabilities", class(probabilities))
  return(probabilities)
}


#' @rdname species_distribution
#' 
#' @export
is_probabilities <- function(x) {
  inherits(x, "probabilities")
}


#  Abundances ----

#' @rdname species_distribution
#'
#' @param round If `TRUE`, the values of `x` are converted to integers.
#' 
#' @export
abundances <- function(
    x,
    round = TRUE,
    names = NULL, 
    weights = NULL, 
    check_arguments = TRUE) {
  
  if (!is.numeric(x)) stop("'x' must be numeric")
  if (any(x < 0)) stop("Species abundances must be positive.")
  if (round) {
    # Add 0.5 before changing mode to round rather than taking the floor
    x <- x + 0.5
    mode(x) <- "integer"
  }
  
  abundances <- species_distribution(
    x,     
    names = names, 
    weights = weights, 
    check_arguments = check_arguments
  )

  class(abundances) <- c("abundances", class(abundances))
  return(abundances)    

}

#' @rdname species_distribution
#' 
#' @export
as_abundances <- function(x, ...) {
  UseMethod("as_abundances")
}


#' @rdname species_distribution
#' 
#' @export
as_abundances.numeric <- function(
    x,
    round = TRUE, 
    ...,
    check_arguments = TRUE) {

  if (check_arguments) check_divent_args()
  if (any(x < 0)) stop("Species abundances must be positive.")
  
  if (round) {
    x <- as.integer(round(x))
  }

  abundances <- as_species_distribution(x, ..., check_arguments = FALSE)

  class(abundances) <- c("abundances", class(abundances))
  return(abundances)
}


#' @rdname species_distribution
#' 
#' @export
as_abundances.matrix <- function(
    x,
    round = TRUE,
    names = NULL, 
    weights = NULL, 
    ...,
    check_arguments = TRUE) {
  
  if (round) {
    # Add 0.5 before changing mode to round rather than taking the floor
    x <- x + 0.5
    mode(x) <- "integer"
  }
  
  abundances <- as_species_distribution.matrix(
    x, 
    names = names, 
    weights = weights,
    ..., 
    check_arguments = check_arguments
  )
  
  class(abundances) <- c("abundances", class(abundances))
  return(abundances)
}


#' @rdname species_distribution
#'
#' @export
as_abundances.data.frame <- function(
    x, 
    ...,
    check_arguments = TRUE) {
  
  abundances <- as_species_distribution.data.frame(
    x, 
    ..., 
    check_arguments = check_arguments
  )
  
  class(abundances) <- c("abundances", class(abundances))
  return(abundances)
}


#' @rdname species_distribution
#' 
#' @export
is_abundances <- function(x) {
  inherits(x, "abundances")
}
