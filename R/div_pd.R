#' Faith's Phylogenetic Diversity of a Community
#'
#' Estimate PD \insertCite{Faith1992}{divent} or FD \insertCite{Petchey2002}{divent} from
#' abundance or probability data and a phylogenetic or functional dendrogram.
#'
#' Estimators to deal with incomplete sampling are not implemented.
#' Use function [div_hill] with argument `q = 0` if they are needed.
#'
#' PD and FD are defined as the total length of the branches of the tree.
#'
#' All species of the `species_distribution` must be found in the tips of the
#' `tree`.
#'
#' @inheritParams check_divent_args
#' @param x An object, that may be a named numeric vector (names are species names)
#' containing abundances or probabilities,
#' or an object of class [abundances] or [probabilities].
#' @param prune What to do when some species are in the tree but not in `x`?
#' If `TRUE`, the tree is pruned to keep species of `x` only.
#' The height of the tree may be changed if a pruned branch is related to the root.
#' If `FALSE` (default), the length of branches of missing species is not summed
#' but the height of the tree is never changed.
#' @param ... Unused.
#'
#' @returns A tibble with the site names, the estimators used and the estimated diversity.
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' # diversity of each community
#' div_pd(paracou_6_abd, tree = paracou_6_taxo)
#'
#' # gamma diversity
#' div_pd(paracou_6_abd, tree = paracou_6_taxo, gamma = TRUE)
#'
#' @name div_pd
NULL


#' @rdname div_pd
#'
#' @export
div_pd <- function(x, tree, ...) {
  UseMethod("div_pd")
}


#' @rdname div_pd
#'
#' @export
div_pd.numeric <- function(
    x,
    tree,
    prune = FALSE,
    as_numeric = FALSE,
    ...,
    check_arguments = TRUE) {

  # Check arguments
  if (any(check_arguments)) {
    check_divent_args()
    if (any(x < 0)) {
      cli::cli_abort("Species probabilities or abundances must be positive.")
    }
  }
  # Check species names
  species_names <- names(x)
  # Prepare the tree
  tree <- as_phylo_divent(tree)
  if (any(check_arguments)) {
    # Check species in the tree
    if (length(setdiff(species_names, rownames(tree$phylo_groups))) != 0) {
      cli::cli_abort("Some species are missing in the tree.")
    }
  }

  # More species in the tree than in x?
  if (prune) {
    species_not_found <- setdiff(tree$phylo$tip.label, names(x))
    if (length(species_not_found) > 0){
      # Prune the tree to keep species in Ps only
      # tree$phylo is the only updated item of tree because others are useless
      tree$phylo <- ape::drop.tip(tree$phylo, species_not_found)
    }
  }

  # Branch lengths
  lengths <- tree$phylo$edge.length
  # Get unnormalized probabilities p(b)
  ltips <- lapply(tree$phylo$edge[, 2], FUN = function(node) tips(tree, node))
  branches <- unlist(
    lapply(ltips, FUN = function(tips.vector) sum(x[tips.vector]))
  )
  # Sum the lengths of branches with abundance > 0
  the_diversity <- sum(lengths[branches != 0])

  # Return
  if (as_numeric) {
    return(the_diversity)
  } else {
    return(
      tibble::tibble_row(
        estimator = "naive",
        order = 0,
        diversity = the_diversity
      )
    )
  }
}


#' @rdname div_pd
#'
#' @export
div_pd.species_distribution <- function(
    x,
    tree,
    prune = FALSE,
    gamma = FALSE,
    as_numeric  = FALSE,
    ...,
    check_arguments = TRUE) {

  # Check arguments
  if (any(check_arguments)) {
    check_divent_args()
    if (any(x < 0)) {
      cli::cli_abort("Species probabilities or abundances must be positive.")
    }
    # Check species names
    col_names <- colnames(x)
    species_names <- col_names[!col_names %in% non_species_columns]
    # Prepare the tree
    tree <- as_phylo_divent(tree)
    # Check species in the tree
    if (length(setdiff(species_names, rownames(tree$phylo_groups))) != 0) {
      cli::cli_abort("Some species are missing in the tree.")
    }
  } else {
    # Prepare the tree
    tree <- as_phylo_divent(tree)
  }

  if (gamma) {
    # Build the metacommunity
    abd <- metacommunity.abundances(
      x,
      as_numeric = TRUE,
      check_arguments = FALSE
    )
    the_diversity <- div_pd.numeric(
      abd,
      tree = tree,
      as_numeric = as_numeric,
      check_arguments = FALSE
    )
    if (!as_numeric) {
      # Add the site column
      the_diversity <- dplyr::bind_cols(
        site = "Metacommunity",
        the_diversity
      )
    }
    return(the_diversity)
   } else {
    # Apply div_pd.numeric() to each site
    div_pd_list <- apply(
      # Eliminate site and weight columns
      x[, !colnames(x) %in% non_species_columns],
      # Apply to each row
      MARGIN = 1,
      FUN = div_pd.numeric,
      # Arguments
      tree = tree,
      prune = prune,
      as_numeric = FALSE,
      check_arguments = FALSE
    )
    # Make a tibble with site, estimator and diversity
    the_diversity <- tibble::tibble(
      # Restore non-species columns
      x[colnames(x) %in% non_species_columns],
      # Coerce the list returned by apply into a dataframe
      do.call(rbind.data.frame, div_pd_list)
    )
    if (as_numeric) {
      return(the_diversity$diversity)
    } else {
      return(the_diversity)
    }
  }
}


#' Names of the tips descending from the node
#'
#' Mimics `geiger::tips()`.
#' Much slower but saves a package dependency.
#'
#' @param node an internal node number of the phylo object.
#'
#' @returns a string: the name of the tip.
#' @noRd
tips <- function(tree, node)
{
  tips_n <- length(tree$phylo$tip.label)
  if (node > tips_n) {
    # internal node numbers start after the last tip
    return(ape::extract.clade(tree$phylo, node)$tip.label)
  } else {
    return(tree$phylo$tip.label[node])
  }
}
