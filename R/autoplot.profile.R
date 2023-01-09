#' Plot Profile Objects
#' 
#' Plot objects of class "profile" produced by [profile_hill] and other 
#' profile functions.
#'
#' @param object An object of class "profile".
#' @param ... Unused.
#' @param main The main title of the plot.
#' @param xlab The label of the x-axis.
#' @param ylab The label of the y-axis.
#' @param shade_color The color of the shaded confidence envelopes.
#' @param alpha The opacity of the confidence envelopes, between 0 (transparent) and 1 (opaque).
#' @param lty The line type of the curves.
#' @param lwd The line width of the curves.
#'
#' @return A [ggplot2::ggplot] object.
#' @export
#'
#' @examples
#' # Diversity profile curve
#' autoplot(profile_hill(mock_3sp_abd))
#' 
autoplot.profile <-  function(
    object, 
    ..., 
    main = NULL,
    xlab = "Order of Diversity",
    ylab = "Diversity",
    shade_color = "grey75",
    alpha = 0.3,
    lty = ggplot2::GeomLine$default_aes$linetype,
    lwd = ggplot2::GeomLine$default_aes$linewidth){
  
  # Add a site column if needed
  if (!"site" %in% colnames(object)) {
    object <- dplyr::mutate(object, site = "Unique site")
  }
  
  # Build the plot
  the_plot <- ggplot2::ggplot(
    object, 
    ggplot2::aes(
      x = .data$order, 
      y = .data$diversity,
      col = .data$site
    )
  )
  
  # Confidence envelope
  if ("sup" %in% colnames(object) & "inf" %in% colnames(object)) {
    the_plot <- the_plot +
      ggplot2::geom_ribbon(
        ggplot2::aes(
          ymin = .data$inf, 
          ymax = .data$sup,
          lty = .data$site
        ), 
        fill = shade_color, 
        alpha = alpha,
        col = shade_color
      )
  }
  
  # Profiles
  the_plot <- the_plot +
    ggplot2::geom_line(linetype = lty, linewidth = lwd) +
    ggplot2::labs(title = main, x = xlab, y = ylab)
  
  return(the_plot)
}
