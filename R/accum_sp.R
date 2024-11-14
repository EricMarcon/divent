#' Spatial Accumulation of Diversity
#'
#' A spatial accumulation is a measure of diversity with respect to the distance from individuals.
#' 
#' Objects of class `accum_sp` contain the value of diversity 
#' (`accum_sp_diversity` objects), entropy (`accum_sp_entropy` objects) or
#' mixing (`accum_sp_mixing` objects) at distances from the individuals.
#' 
#' These objects are lists: 
#' 
#' - `X` contains the [dbmss::wmppp] point pattern,
#' - `accumulation` is a 3-dimensional array, with orders of diveristy in rows,
#' neighborhood size (number of points or distance) in columns and a single slice
#' for the observed entropy, diversity or mixing.
#' - `neighborhoods` is a similar 3-dimensional array with one slice per point
#' of `X`.
#' 
#' They can be plotted or mapped.
#' @aliases accum_sp accum_sp_entropy accum_sp_diversity accum_sp_mixing
#' 
#' @inheritParams check_divent_args
#' 
#' @name accum_sp
NULL

#' @rdname accum_sp
#' 
#' @param x an `accum_sp` object.
#' @param ... Additional arguments to be passed to [plot], or, in `plot_map()`,
#' to [spatstat.explore::bw.smoothppp] and [spatstat.explore::density.ppp] to 
#' control the kernel smoothing and to [spatstat.geom::plot.im] to plot the image.
#' @param type plotting parameter. Default is "l".
#' @param main main title of the plot.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param ylim limits of the Y-axis, as a vector of two numeric values.
#' @param show_h0 if `TRUE`, the values of the null hypothesis are plotted.
#' @param line_width width of the Diversity Accumulation Curve line.
#' @param col_shade The color of the shaded confidence envelope.
#' @param col_border The color of the borders of the confidence envelope.
#' 
#' @returns `plot.accum_sp()` returns `NULL`.
#' @export
#'
plot.accum_sp <- function(
    x, 
    ..., 
    q = dimnames(x$accumulation)$q[1], 
    type = "l", 
    main = "accumulation of ...", 
    xlab = "Sample size...", 
    ylab = "Diversity...", 
    ylim = NULL,
    show_h0 = TRUE, 
    line_width = 2, 
    col_shade = "grey75", 
    col_border = "red")  {
  
  # Prepare the parameters
  h <- accum_sp_plot_helper(x, q, main, xlab, ylab, ylim)
  
  # Prepare the plot
  plot(
    x = dimnames(x$accumulation)[[2]], 
    y = as.numeric(x$accumulation[h$q_row, , 1]), 
    ylim = c(h$ymin, h$ymax),
    type = h$type, 
    main = h$main, 
    xlab = h$xlab, 
    ylab = h$ylab
  )
  
  if (dim(x$accumulation)[3] == 4) {
    # Confidence envelope is available
    graphics::polygon(
      x = c(rev(dimnames(x$accumulation)[[2]]), dimnames(x$accumulation)[[2]]), 
      y = c(rev(x$accumulation[h$q_row, , 4]), x$accumulation[h$q_row, , 3]), 
      col = col_shade, 
      border = FALSE
    )
    # Add red lines on borders of polygon
    graphics::lines(
      x = dimnames(x$accumulation)[[2]], 
      y = x$accumulation[h$q_row, , 4], 
      col = col_border, 
      lty = 2
    )
    graphics::lines(
      x = dimnames(x$accumulation)[[2]], 
      y = x$accumulation[h$q_row, , 3], 
      col = col_border, 
      lty = 2
    )
    # Redraw the SAC
    graphics::lines(
      x = dimnames(x$accumulation)[[2]], 
      y = x$accumulation[h$q_row, , 1], 
      lwd = line_width, 
      ...
    )
    
    # H0
    if (show_h0) {
      graphics::lines(
        x = dimnames(x$accumulation)[[2]], 
        y = x$accumulation[h$q_row, , 2], 
        lty = 2
      )      
    } 
  }
  
  return(invisible(NULL))
}


#' @rdname accum_sp
#'
#' @param object an `accum_sp` object.
#'
#' @returns `autoplot.accum_sp()` returns a [ggplot2::ggplot] object.
#' @export
#'
autoplot.accum_sp <- function(
    object, 
    ..., 
    q = dimnames(object$accumulation)$q[1],
    main = "Accumulation of ...", 
    xlab = "Sample size...", 
    ylab = "Diversity...", 
    ylim = NULL, 
    show_h0 = TRUE, 
    col_shade = "grey75", 
    col_border = "red")   {
  
  # Prepare the parameters
  h <- accum_sp_plot_helper(object, q, main, xlab, ylab, ylim)
  
  # Prepare the data
  df <- data.frame(
    x = as.numeric(dimnames(object$accumulation)[[2]]), 
    y = object$accumulation[h$q_row, , 1]
  )
  if (dim(object$accumulation)[3] == 4) {
    # Confidence envelope is available
    df$low <- object$accumulation[h$q_row, , 3]
    df$high <- object$accumulation[h$q_row, , 4]
    if (show_h0) df$H0 <- object$accumulation[h$q_row, , 2]
  }
  
  # Prepare the plot
  the_plot <- ggplot2::ggplot(
    data = df, 
    ggplot2::aes(x = .data$x, y = .data$y)
  ) +
    ggplot2::geom_line() +
    ggplot2::labs(title = h$main, x = h$xlab, y = h$ylab)
  
  if (dim(object$accumulation)[3] == 4) {
    the_plot <- the_plot +
      ggplot2::geom_ribbon(
        ggplot2::aes(
          ymin = .data$low, 
          ymax = .data$high
        ), 
        fill = col_shade, 
        alpha = 0.5) +
      # Add red lines on borders of polygon
      ggplot2::geom_line(
        ggplot2::aes(y = .data$low), 
        colour = col_border, 
        linetype = 2
      ) +
      ggplot2::geom_line(
        ggplot2::aes(y = .data$high), 
        colour = col_border, 
        linetype = 2
      )
    
    # H0
    if (show_h0) {
      the_plot <- the_plot +
        ggplot2::geom_line(ggplot2::aes(y = .data$H0), linetype = 2)
    }
  }
  return(the_plot)
}


#' @rdname accum_sp
#' 
#' @param accum an object to map.
#' @param neighborhood The neighborhood size, i.e. the number of neighbors or the distance to consider.
#' @param sigma the smoothing bandwidth. 
#' The standard deviation of the isotropic smoothing kernel. 
#' Either a numerical value, or a function that computes an appropriate value of sigma.
#' @param allow_jitter if `TRUE`, duplicated points are jittered to avoid their
#' elimination by the smoothing procedure.
#' @param weighted if `TRUE`, the weight of the points is used by the smoothing
#' procedure.
#' @param adjust force the automatically selected bandwidth to be multiplied 
#' by `adjust`. 
#' Setting it to values lower than one (1/2 for example) will sharpen the estimation.
#' @param dim_x the number of columns (pixels) of the resulting map, 128 by default.
#' @param dim_y the number of rows (pixels) of the resulting map, 128 by default.
#' @param col the colors of the map. See [spatstat.geom::plot.im] for details.
#' @param contour if `TRUE`, contours are added to the map.
#' @param contour_levels the number of levels of contours.
#' @param contour_col the color of the contour lines.
#' @param points if `TRUE`, the points that brought the data are added to the map.
#' @param pch the symbol used to represent points.
#' @param point_col The color of the points.
#' Standard base graphic arguments such as `main` can be used.
#' 
#' @returns `plot_map` returns a [spatstat.geom::im] object that can be used to produce 
#' alternative maps.
# 
#' @export
#' 
#' @examples
#' # Generate a random community
#' X <- rspcommunity(1, size = 50, species_number = 10)
#' # Calculate the species accumulation curve
#' accum <- accum_sp_hill(X, orders = 0, r = c(0, 0.2), individual = TRUE)
#' # Plot the local richness at distance = 0.2
#' plot_map(accum, q = 0, neighborhood = 0.2)
#' 
plot_map <- function(
    accum, 
    q = dimnames(accum$accumulation)$q[1], 
    neighborhood = dplyr::last(colnames(accum$neighborhoods)), 
    sigma = spatstat.explore::bw.scott(accum$X, isotropic = TRUE), 
    allow_jitter = TRUE,
    weighted = FALSE, 
    adjust = 1, 
    dim_x = 128, 
    dim_y = 128, 
    col = grDevices::terrain.colors(256), 
    contour = TRUE,
    contour_levels = 10, 
    contour_col = "dark red",
    points = FALSE, 
    pch = 20, 
    point_col = "black",
    ..., 
    check_arguments = TRUE) {
  
  if (any(check_arguments)) {
    check_divent_args()
  }
  if (is.null(dim(accum$neighborhoods))) {
    stop("The Accumulation object does not contain individual data to plot.")
  }
    
  # Jitter
  if (allow_jitter) {
    # Find duplicates
    the_dups <- spatstat.geom::duplicated.ppp(accum$X, rule = "unmark")
    if (sum(the_dups) > 0) {
      # Extract the duplicates and jitter them
      the_dups.wmppp <- spatstat.geom::rjitter(accum$X[the_dups])
      # Put the coordinates back into the original wmppp
      accum$X$x[the_dups] <- the_dups.wmppp$x
      accum$X$y[the_dups] <- the_dups.wmppp$y
    }
  }
  
  # Convert numeric values of q and Neighborhood into their index
  q_row <- which(as.numeric(rownames(accum$neighborhoods)) == q)
  nbd_col <- which(as.numeric(colnames(accum$neighborhoods)) == neighborhood)
  # Verify that values exist: if which() did not match, we get integer(0) for q or neighborhood
  # then data is of length 0.
  if (length(accum$neighborhoods[q_row, , ]) == 0) {
    stop("Incorrect q.")
  }
  if (length(accum$neighborhoods[, nbd_col, ]) == 0) {
    stop("Incorrect neighborhood.") 
  }
  
  # Detect points with NA values
  is_not_na <- !is.na(accum$neighborhoods[q_row, nbd_col, ])
  
  # Weights
  if (weighted) {
    the_weights <- spatstat.geom::marks(accum$X)[is_not_na]
  } else {
    the_weights <- rep(1, spatstat.geom::npoints(accum$X))
  }
  
  # Prepare the ppp to plot
  the_ppp <- spatstat.geom::unmark(accum$X)
  spatstat.geom::marks(the_ppp) <- accum$neighborhoods[q_row, nbd_col, ]
  the_ppp <- the_ppp[is_not_na]
  class(the_ppp) <- "ppp"
  
  # Image
  the_image <- spatstat.explore::Smooth.ppp(
    the_ppp, 
    sigma = sigma, 
    ..., 
    weights = the_weights, 
    adjust = adjust,
    dimyx = c(dim_y, dim_x)
  )
  
  plot(the_image, col = col, ...)
  if (contour) {
    graphics::contour(
      the_image, 
      add = TRUE, 
      nlevels = contour_levels, 
      col = contour_col
    )
  }
  if(points) {
    graphics::points(
      x = the_ppp$x, 
      y = the_ppp$y, 
      pch = pch, 
      col = point_col
    )
  }

  # Return the image for further processing
  return(invisible(the_image))
}
