#' Class phylo_divent
#' 
#' Methods for dendrograms of class "phylo_divent".
#' 
#' `as_phylo_divent` calculates cuts and intervals of a phylogenetic tree and makes
#' it available both in [stats::hclust] and [ape::phylo] formats.
#' The conversion preprocesses the tree: it calculates cuts so that the tree
#' can be reused efficiently by phylodiversity functions.
#'
#' @param x An object of class "phylo_divent".
#' @param ... Arguments passed to [plot.phylo].
#'
#' @return `as_phylo_divent` returns a phylogenetic tree that is an object of 
#' class "phylo_divent".
#'
#' @examples
#' # Paracou plot 6 species taxonomy 
#' plot(paracou_6_taxo, type="fan", show.tip.label=FALSE)
#' 
#' @name phylo_divent
NULL


#' @rdname phylo_divent
#' 
#' @param tree An object of class [ape::phylo], [ade4::phylog] or [stats::hclust].
#'
#' @export
as_phylo_divent <- function(tree) {
  # The tree may be NULL or already processed
  if (is.null(tree) | inherits(tree, "phylo_divent")) return (tree)
  
  # tree must be either a phylog, phylo or a hclust object
  if (inherits(tree, "phylog")) {
    # Build an hclust object to use cutree later. 
    # Distances in $Wdist are actually sqrt(2*distance)
    # Caution: Distances in hclust count full branch lengths between species, 
    # i.e. twice the ultrametric distance. See ?as.phylo.hclust
    tree.hclust <- stats::hclust(tree$Wdist^2 / 2, "average")
    # build a phylo object
    tree.phylo <- ape::as.phylo.hclust(tree.hclust)
    # Double edge.lengths to correct as.phylo.hclust
    tree.phylo$edge.length <- 2 * tree.phylo$edge.length
  } else {
    if (inherits(tree, "phylo")) {
      tree.phylo <- tree
      # Build an hclust object to use cutree later.
      # Edge lengths are multiplied by 2 during the conversion. 
      # Divide by 2 before that.
      tree$edge.length <- tree$edge.length / 2
      tree.hclust <- ape::as.hclust.phylo(tree)
    } else {
      if (inherits(tree, "hclust")) {
        # Caution: Distances in hclust count full branch lengths between species,
        # i.e. twice the ultramtetric distance
        tree.hclust <- tree
        # build a phylo object to use $droot later
        tree.phylo <- ape::as.phylo.hclust(tree)
        # Double edge.lengths to correct as.phylo.hclust
        tree.phylo$edge.length <- 2 * tree.phylo$edge.length
      } else {
        stop("tree must be an object of class phylo, phylog or hclust")
      }
    }
  }
  
  # Calculate distances between nodes and leaves
  dist_from_leaves <- ape::branching.times(tree.phylo)
  # Get a sorted list of cuts (eliminate leaves)
  cuts <- sort(unique(dist_from_leaves))
  # Calculate intervals between cuts (add 0 to cuts to get the first interval)
  intervals <- diff(c(0, cuts))
  # Eliminate 0 intervals (happen when a node contains more than 2 tips), 
  # including rounding errors
  rounding_error <- max(dist_from_leaves) * 10 * .Machine$double.eps
  cuts <- cuts[intervals > rounding_error]
  intervals <- intervals[intervals > rounding_error]
  
  # Format and return
  the_tree <- list(
    phylo     = tree.phylo,
    hclust    = tree.hclust,
    height    = cuts[length(cuts)],
    cuts      = cuts,
    intervals = intervals
  )
  class(the_tree) <- "phylo_divent"
  return(the_tree)
}


#' @rdname phylo_divent
#'
#' @export
is_phylo_divent <- function (x) {
  inherits(x, "phylo_divent")
}


#' @rdname phylo_divent
#'
#' @importFrom graphics plot
#' @export
plot.phylo_divent <- function (x, ...) {  
  # Plot the phylo object
  if (requireNamespace("ape")) {
    plot(x$phylo, ...)
  } else {
    warning("The package 'ape' is required to plot the tree.")
  }
}
