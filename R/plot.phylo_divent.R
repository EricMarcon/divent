#' Plot phylo_divent Objects
#' 
#' Plot objects of class "phylo_divent" produced by [as_phylo_divent], that are
#' phylogenetic trees.
#'
#' @param x 
#' @param ... Arguments passed to [stats::plot.dendrogram].
#'
#' @importFrom graphics plot
#' @export
#'
#' @examples
#' # Paracou plot 6 species taxonomy
#' tree <- as_phylo_divent(paracou_6_taxo)
#' plot(tree, type="fan", show.tip.label=FALSE)
#' 
plot.phylo_divent <- function (x, ...) {  
  # Plot the hclust object as a dendrogram
  plot(stats::as.dendrogram(x$hclust), ...)
}
