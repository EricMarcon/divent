#' Plot Profile Objects
#' 
#' Plot objects of class "species_distribution" produced by 
#' [species_distribution] and similar functions.
#' 
#' @param x An object.
#' @param ... Additional arguments to be passed to [plot]. Unused elsewhere.
#'
#' @name plot.species_distribution
NULL

#' @rdname plot.species_distribution
#'
#' @param type The type of plot. "RAC" (Rank-abundance curve, or Whittaker plot)
#' or "Metacommunity" to represent species abundances of each community along 
#' with those of the metacommunity.
#' @param fit_rac If `TRUE`, estimate a theoretical distribution and fit the data with it.
#' RAC plot only.
#' @param distribution The distribution of species abundances.
#' May be "lnorm" (log-normal), "lseries" (log-series), "geom" (geometric) or 
#' "bstick" (broken stick).
#' RAC plot only.
#' @param ylog If `TRUE`, the Y-axis is in log-scale.
#' RAC plot only.
#' @param main The title of the plot.
#' @param xlab The label of the X-axis.
#' RAC plot only.
#' @param ylab The  label of the Y-axis.
#' @param palette The name of a color palette, recognized by [RColorBrewer::brewer.pal].
#' RAC plot only.
#'
#' @importFrom base plot
#' @export
plot.species_distribution <- function(
    x, 
    type = c("RAC", "Metacommunity"),
    ..., 
    fit_rac = FALSE,
    distribution = c("lnorm", "lseries", "geom", "bstick"),
    ylog = "y", 
    main = NULL, 
    xlab = "Rank", 
    ylab = NULL,
    palette = "Set1") {
  
  type <- match.arg(type)
  
  if (type == "RAC") {
    # Whittaker plot ----
    
    # Prepare ylab
    if (is.null(ylab)) {
      if (is_probabilities(x)) {
        ylab <- "Probability"
      } else {
        ylab <- "Abundance" 
      }
    }
    
    # Find the max number of species in all communities
    abundances <- x[, !colnames(x) %in% c("site", "weight")] 
    s_obs_max <- max(rowSums(abundances > 0))
    
    # Prepare the plot: X and Y ranges
    base::plot(
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
      abd <- x[community, !colnames(x) %in% c("site", "weight")]
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
  } else if (type == "Metacommunity") {
    # Metacommunity plot ----
    
    # Prepare ylab
    if (is.null(ylab)) {
      ylab <- "Species frequencies"
    }
    # Prepare data: community probabilities
    x.probabilities <- probabilities.abundances(
      x, 
      estimator = "naive", 
      check_arguments = FALSE
    )
    prob_communities <- t(
      x.probabilities[, !colnames(x.probabilities) %in% non_species_columns]
    )
    # Normalize weights (that are the widths of community bars)
    weights <- x$weight / sum(x$weight)
    # Metacommunity probabilities
    x.metacommunity <- metacommunity(x)
    abd_metacommunity <- as.numeric(
      x.metacommunity[1, !colnames(x.metacommunity) %in% non_species_columns]
    )
    prob_metacommunity <- abd_metacommunity / sum(abd_metacommunity)
    
    # Plot
    graphics::barplot(
      cbind(
        prob_communities,
        rep(0, length(abd_metacommunity)),
        prob_metacommunity 
      ),
      beside = FALSE,
      width = c(weights, .5, 1),
      names.arg = c(x.probabilities$site, "", "Metacommunity"),
      main = main,
      ylab = ylab,
      ...
    )
  }
  
}
#' @export
base::plot


#' @rdname plot.species_distribution
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
      object[, !colnames(object) %in% non_species_columns] > 0
    )
  )
  
  # Prepare the plot
  the_plot <- ggplot2::ggplot() +
    # Make X axis start at 1
    ggplot2::scale_x_continuous(labels = rlang::as_function(~ .x + 1)) +
    ggplot2::labs(title = main, x = xlab, y = ylab)
  
  # Prepare the data to plot
  the_data <- data.frame(site = NULL, rank = NULL, abd = NULL)
  if (fit_rac) {
    the_model <- data.frame(site = NULL, rank = NULL, abd = NULL)
  }
  
  # Loop in communities to build the plot
  for (community in seq_len(nrow(object))) {
    # Extract the abundances of the community
    abd <- object[community, !colnames(object) %in% non_species_columns]
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
