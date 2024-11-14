#' Plot a map
#'
#' Map Spatialized Communities
#' 
#' Maps of interpolated values from points of an [accum] produced e.g. by the function [DivAccum] object are produced.
#' 
#' @inheritParams check_divent_args
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
#' @param ... further arguments passed to [spatstat.explore::bw.smoothppp] and 
#' [spatstat.explore::density.ppp] to control the kernel smoothing and 
#' to [spatstat.geom::plot.im] to plot the image.
#' Standard base graphic arguments such as `main` can be used.
#' 
#' @returns A [spatstat.geom::im] object that can be used to produce 
#' alternative maps.
# 
#' @export
#' 
#' @examples
#' # Generate a random community
#' X <- rspcommunity(1, size = 50, species_number = 10)
#' # Calculate the species accumulation curve
#' accum <- accum_sp_hill(X, orders = 0, r = c(0, 0.2), individual=TRUE)
#' # Plot the local richness at distance = 0.2
#' plot_map(accum, q = 0, neighborhood = 0.2)
#' 
plot_map <- function(
    accum, 
    q = dimnames(accum$Accumulation)$q[1], 
    neighborhood = dplyr::last(colnames(accum$Neighborhoods)), 
    sigma = spatstat.explore::bw.scott(accum$X, isotropic = TRUE), 
    allow_jitter = TRUE,
    weighted = FALSE, 
    adjust = 1, 
    dim_x = 128, 
    dim_y = 128, 
    col = terrain.colors(256), 
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
  if (is.null(dim(accum$Neighborhoods))) {
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
  q_row <- which(as.numeric(rownames(accum$Neighborhoods)) == q)
  neighborhood <- which(as.numeric(colnames(accum$Neighborhoods)) == neighborhood)
  # Verify that values exist: if which() did not match, we get integer(0) for q or neighborhood
  # then data is of length 0.
  if (length(accum$Neighborhoods[q, , ]) == 0) {
    stop("Incorrect q.")
  }
  if (length(accum$Neighborhoods[, neighborhood, ]) == 0) {
    stop("Incorrect neighborhood.") 
  }
  
  # Detect points with NA values
  is_not_na <- !is.na(accum$Neighborhoods[q, neighborhood, ])
  
  # Weights
  if (weighted()) {
    the_weights <- spatstat.geom::marks(accum$X)[!is_na]
  } else {
    the_weights <- rep(1, spatstat.geom::npoints(accum$X))
  }
  
  # Prepare the ppp to plot
  the_ppp <- spatstat.geom::unmark(accum$X)
  spatstat.geom::marks(the_ppp) <- accum$Neighborhoods[q, neighborhood, ]
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
