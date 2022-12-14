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
#' `species_distribution` objects objects can be plotted by [plot] and [autoplot].
#'
#' @inheritParams check_divent_args
#' @param x An object.
#' @param ... Additional arguments to be passed to [plot]. Unused elsewhere.
#' 
#' @references
#' \insertAllCited{}
#' 
#' @examples
#' # Paracou data is a tibble
#' paracou_6_abd
#' # Class
#' class(paracou_6_abd)
#' is_species_distribution(paracou_6_abd)
#' # Whittaker plot fitted by a log-normal distribution
#' autoplot(paracou_6_abd[1,], fit_rac = TRUE, distribution = "lnorm")
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
  if (any(check_arguments)) check_divent_args()
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
    the_distribution <- tibble::tibble_row(
      site = names,
      weight = sum(x),
      as.data.frame(t(x))
    )
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
    the_distribution <- tibble::as_tibble(x, rownames = "site")
    ### Rows: site names = names or matrix row names or default ----
    if (!is.null(names)) {
      # site = names if the size matches
      if (length(names) == nrow(x)) {
        the_distribution$site <- names
      } else {
        stop("The length of 'names' must match the number of lines of the data matrix.")
      }
    } else {
      # names is null...
      if (is.null(row.names(x))) {
        # ...and no row names: set default names such as site_1
        the_distribution$site <- paste(
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
        the_distribution <- tibble::add_column(
          the_distribution, 
          weight = rowSums(x),
          .after = "site"
        )
      } else {
        stop("The length of 'weights' must match the number of lines of the data matrix.")
      }
    } else {
      # Weights are the number of individuals
      the_distribution <- tibble::add_column(
        the_distribution, 
        weight = rowSums(x),
        .after = "site"
      )
    }
  }
  
  # Set the class and return ----
  class(the_distribution) <- c("species_distribution", class(the_distribution))
  return(the_distribution)
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
  
  if (any(check_arguments)) check_divent_args()
  # Check the data
  if (any(x < 0)) stop("All numeric values of the dataframe must be positive.")
  
  # Build a tibble
  the_distribution <- tibble::as_tibble(x)
  
  # The first column should be "site"
  if (!"site" %in% colnames(the_distribution)) {
    the_distribution <- tibble::add_column(
      the_distribution, 
      site = paste(
        "site", 
        formatC(
          seq_len(nrow(the_distribution)), 
          width = ceiling(log10(nrow(the_distribution))), 
          flag = "0"
        ),
        sep = "_"
      ),
      .before = 1
    )
  }
  
  # The second column should be "weight"
  if (!"weight" %in% colnames(the_distribution)) {
    the_distribution <- tibble::add_column(
      the_distribution, 
      weight = rowSums(the_distribution[colnames(the_distribution) != "site"]),
      .after = "site"
    )
  }

  # Set the class and return
  class(the_distribution) <- c("species_distribution", class(the_distribution))
  return(the_distribution)
}


#' @rdname species_distribution
#'
#' @export
is_species_distribution <- function(x) {
  inherits(x, "species_distribution")
}


#' @rdname species_distribution
#'
#' @param fit_rac If `TRUE`, estimate a theoretical distribution and fit the data with it.
#' @param distribution The distribution of species abundances.
#' May be "lnorm" (log-normal), "lseries" (log-series), "geom" (geometric) or 
#' "bstick" (broken stick).
#' @param ylog If `TRUE`, the Y-axis is in log-scale.
#' @param main The title of the plot.
#' @param xlab The label of the X-axis.
#' @param ylab The  label of the Y-axis.
#' @param palette The name of a color palette, recognized by [RColorBrewer::brewer.pal].
#'
#' @importFrom graphics plot
#' @export
plot.species_distribution <- function(
    x, 
    ..., 
    fit_rac = FALSE,
    distribution = c("lnorm", "lseries", "geom", "bstick"),
    ylog = "y", 
    main = NULL, 
    xlab = "Rank", 
    ylab = NULL,
    palette = "Set1") {
  
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
    log = ylog, 
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
  
  # Color palette, min number of colors is 3
  cols <- RColorBrewer::brewer.pal(max(nrow(x), 3), name = palette)
  
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
    if (fit_rac) {
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
    graphics::legend(
      "topright",
      inset = .02,
      legend = x$site,
      col = cols,
      lty = 1,
      pch = 1
    )
  }
}
#' @export
graphics::plot


#' @rdname species_distribution
#'
#' @param object An object of class [species_distribution].
#' @param pch The plotting characters. See [graphics::points].
#' @param cex The character expansion (size) of the points. See [graphics::points].
#'
#' @importFrom ggplot2 autoplot
#' @importFrom rlang .data
#' @export
autoplot.species_distribution <- function(
    object, 
    ..., 
    fit_rac = FALSE,
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
  if (fit_rac) {
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
    if (fit_rac) {
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
  if (fit_rac) {
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
#' @export
ggplot2::autoplot


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
  the_probabilities <- as_species_distribution(
    prob,
    check_arguments = check_arguments
  )

  class(the_probabilities) <- c("probabilities", class(the_probabilities))
  return(the_probabilities)
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
  
  # Calculate probabilities by row
  prob <- x / rowSums(x)
  
  # Build the species distribution
  the_probabilities <- as_species_distribution.matrix(
    prob, 
    names = names, 
    weights = weights,
    check_arguments = check_arguments
  )
  
  # Set the class
  class(the_probabilities) <- c("probabilities", class(the_probabilities))
  return(the_probabilities)
}


#' @rdname species_distribution
#'
#' @export
as_probabilities.data.frame <- function(
    x, 
    ...,
    check_arguments = TRUE) {
  
  # Build the species distribution to add site and weight columns if needed
  abundances <- as_species_distribution.data.frame(
    x, 
    check_arguments = check_arguments
  )
  
  # Select species columns
  species_columns <- !(colnames(abundances) %in% c("site", "weight"))
  # Extract abundances
  abd <- as.matrix(abundances[, species_columns])
  # Normalize them
  prob <- abd / rowSums(abd)
  # Build the tibble
  the_probabilities <- cbind(
    data.frame(site = abundances$site, weight = abundances$weight),
    as.data.frame(prob)
  )
  # Restore exact species names (spaces may have been transformed into "_")
  colnames(the_probabilities[, species_columns]) <- colnames(abundances[, species_columns])
  
  # Build the species distribution again for the classes
  the_probabilities <- as_species_distribution.data.frame(
    the_probabilities,
    check_arguments = FALSE
  )
  
  # Set the class
  class(the_probabilities) <- c("probabilities", class(the_probabilities))
  return(the_probabilities)
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

  if (any(check_arguments)) check_divent_args()
  if (any(x < 0)) stop("Species abundances must be positive.")
  
  if (round) {
    x <- as.integer(round(x))
  }

  the_abundances <- as_species_distribution(x, check_arguments = FALSE)

  class(the_abundances) <- c("abundances", class(the_abundances))
  return(the_abundances)
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
  
  the_abundances <- as_species_distribution.matrix(
    x,
    names = names,
    weights = weights,
    check_arguments = check_arguments
  )
  
  class(the_abundances) <- c("abundances", class(the_abundances))
  return(the_abundances)
}


#' @rdname species_distribution
#'
#' @export
as_abundances.data.frame <- function(
    x, 
    ...,
    check_arguments = TRUE) {
  
  the_abundances <- as_species_distribution.data.frame(
    x,
    check_arguments = check_arguments
  )
  
  class(the_abundances) <- c("abundances", class(the_abundances))
  return(the_abundances)
}


#' @rdname species_distribution
#' 
#' @export
is_abundances <- function(x) {
  inherits(x, "abundances")
}
