#' Plot phylo_divent Objects
#'
#' Plot objects of class "phylo_divent" produced by [as_phylo_divent], that are
#' phylogenetic trees.
#'
#' @param x An object of class "phylo_divent".
#' @param ... Arguments passed to [stats::plot.dendrogram].
#'
#' @returns `NULL`. Called for side effects.
#'
#' @export
#'
#' @examples
#' # Paracou plot 6 species taxonomy
#' tree <- as_phylo_divent(paracou_6_taxo)
#' plot(tree, leaflab = "none")
#'
plot.phylo_divent <- function(x, ...) {
  # Plot the hclust object as a dendrogram
  plot(stats::as.dendrogram(x$hclust), ...)
}
