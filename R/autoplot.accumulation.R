#' Plot Accumulation Objects
#' 
#' Plot objects of class "accumulation" produced by [accum_hill] and other 
#' accumulation functions.
#'
#' @param object An object of class "accumulation".
#' @param ... 
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
#' # Species accumulation curve
#' autoplot(div_accum(mock_3sp_abd))
#' 
autoplot.accumulation <-  function(
    object, 
    ..., 
    main = NULL,
    xlab = "Sample Size",
    ylab = NULL,
    shade_color = "grey75",
    alpha = 0.3,
    lty = ggplot2::GeomLine$default_aes$linetype,
    lwd = ggplot2::GeomLine$default_aes$linewidth){
  
  # Add a site column if needed
  if (!"site" %in% colnames(object)) {
    object <- dplyr::mutate(object, site = "Unique site")
  }
  
  # Build the plot
  if ("diversity" %in% colnames(object)) {
    the_plot <- ggplot2::ggplot(
      object, 
      ggplot2::aes(
        x = .data$level, 
        y = .data$diversity,
        col = .data$site
      )
    )
  } else {
    the_plot <- ggplot2::ggplot(
      object, 
      ggplot2::aes(
        x = .data$level, 
        y = .data$entropy,
        color = .data$site
      )
    )
  }
  
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
  
  # y-axis label
  if (is.null(ylab)) {
    if ("diversity" %in% colnames(object)) {
      ylab <- "Diversity"
    } else {
      ylab <- "Entropy"
    }
  }
  
  # Accumulations
  the_plot <- the_plot +
    ggplot2::geom_line(linetype = lty, linewidth = lwd) +
    ggplot2::labs(title = main, x = xlab, y = ylab)
  
  # Actual value
  actual <- dplyr::filter(object, .data$estimator == "Sample")
  if (nrow(actual) > 0) {
    if ("diversity" %in% colnames(object)) {
      the_plot <- the_plot +
        ggplot2::geom_hline(
          data = actual, 
          mapping = ggplot2::aes(yintercept = .data$diversity, color = .data$site), 
          linetype = 2
        )
    } else {
      the_plot <- the_plot +
        ggplot2::geom_hline(
          data = actual, 
          mapping = ggplot2::aes(yintercept = .data$entropy, color = .data$site), 
          linetype = 2
        )
    }
    the_plot <- the_plot +
      ggplot2::geom_vline(
        data = actual, 
        mapping = ggplot2::aes(xintercept = .data$level, color = .data$site), 
        linetype = 2
      ) 
  }
  
  return(the_plot)
}
