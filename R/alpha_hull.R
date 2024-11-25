#' alpha-shape calculation
#'
#' Calculate a window containing all points of a point pattern.
#' The window is not convex but as close as possible to the points.
#'
#' The typical use of this function is to define a narrow window around
#' a point pattern that has been created with a default, rectangle window.
#'
#' The window is built by the [alphahull::ashape()] function and then transformed
#' into a [spatstat.geom::owin.object].
#' The `alpha` parameter determines the smallest size of zones excluded from the window.
#' If it is not specified, a first attempt is 1/256 of the diameter of the existing window of `X`.
#' If the shape cannot be calculated, `alpha` is doubled and a new attempt is made.
#'
#' @param X a planar point pattern ([spatstat.geom::ppp.object]).
#' @param alpha a smoothing parameter to delimit concave polygons.
#'
#' @return A window, i.e. a [spatstat.geom::owin.object].
#' @export
#'
#' @seealso [spatstat.geom::convexhull]
#'
#' @examples
#' # Simulate a point pattern
#' if (require(spatstat.random)) {
#'   X <- rpoispp(50)
#'   plot(X)
#'   # Calculate its border
#'   X$window <- alphahull(X)
#'   plot(X)
#' }
alphahull <- function(X, alpha = NULL) {

  if (!inherits(X, "ppp")) {
    stop("X must be a plana point pattern (ppp object)")
  }
  if (!is.null(alpha)) {
    if (!is.numeric(alpha)) {
      stop("alpha must be numeric")
    }
    if (length(alpha) > 1) {
      stop("alpha must be a single number")
    }
    if (alpha <= 0) {
      stop("alpha must be a positive number.")
    }
  }
  # At least 3 points are needed
  if (X$n < 3)
    return(X$window)
  if (X$n == 3)
    # Window is a convex hull
    return(spatstat.geom::convexhull.xy(X$x, X$y))

  # Prepare
  is_validated_alpha <- FALSE
  # Unique points are needed
  xy <- unique(data.frame(x = X$x, y = X$y))
  # Size of the window
  the_diameter <- spatstat.geom::diameter(spatstat.geom::Window(X))
  # Use alphahull::ashape to obtain a concave hull.
  # Parameter alpha must be small for accuracy (argument alpha)
  if (is.null(alpha)) {
    alpha <- the_diameter / 256
  }
  # but large enough to be able to build a correct graph from the points:
  # alpha is multiplied by 2 until success.
  while (!is_validated_alpha) {
    alpha_shape <- alphahull::ashape(xy, alpha = alpha)
    # Convert alpha shape into polygon (https://rpubs.com/geospacedman/alphasimple)
    # Make a graph with edges, library igraph
    graph_alpha_shape <- igraph::graph.edgelist(
      cbind(
        as.character(alpha_shape$edges[, "ind1"]),
        as.character(alpha_shape$edges[, "ind2"])),
        directed = FALSE
      )
    # Tests: the graph must be connected and circular. If it is not, increase alpha.
    the_error <- ""
    if (alpha_shape$length == 0) {
      the_error <- "No edges in alpha shape"
    } else if (!igraph::is.connected(graph_alpha_shape)) {
      the_error <- "Graph not connected"
    } else if (any(igraph::degree(graph_alpha_shape) != 2)) {
      the_error <- "Graph not circular"
    } else if (igraph::clusters(graph_alpha_shape)$no > 1) {
      the_error <- "Graph composed of more than one circle"
    }
    if (the_error == "") {
      is_validated_alpha <- TRUE
    } else {
      if (alpha > the_diameter) {
        # Unable to make a circular graph: give up.
        warning(paste("Unable to build an alpha hull:", the_error))
      }
      else # Try to double alpha
        alpha <- 2 * alpha
    }
  }

  # Eliminate the first node to destroy circularity
  graph_cut <- graph_alpha_shape - igraph::E(graph_alpha_shape)[1]
  # Find chain end points
  ends <- names(which(igraph::degree(graph_cut) == 1))
  path <- igraph::get.shortest.paths(graph_cut, ends[1], ends[2])[[1]]
  # This is an index into the points
  X_path <- as.numeric(igraph::V(graph_cut)[unlist(path)]$name)
  # Join the ends to restore circularity
  X_path = c(X_path, X_path[1])

  # Get the points from the ashape object, make an owin.
  # Manage reverse by tryCatch
  the_window <- tryCatch(
    spatstat.geom::owin(
      poly = list(
        x = alpha_shape$x[X_path, ][, 1],
        y = alpha_shape$x[X_path, ][, 2]
      )
    ),
    error = function(e) {
      # Error if the polygon is traversed clockwise
      spatstat.geom::owin(
        poly = list(
          x = alpha_shape$x[rev(X_path), ][, 1],
          y = alpha_shape$x[rev(X_path), ][, 2]
        )
      )
    }
  )
  return(the_window)
}
