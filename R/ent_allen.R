#' Allen et al.'s Phylogenetic Entropy of a Community
#'
#' Estimate entropy \insertCite{Allen2009}{divent} from
#' abundance or probability data and a phylogenetic or functional dendrogram.
#'
#' Estimators to deal with incomplete sampling are not implemented.
#' Use function [ent_phylo] with argument if they are needed.
#'
#' The phylogenetic entropy is calculated following
#' \insertCite{Allen2009:textual}{divent} for order \eqn{q=1} and
#' \insertCite{Leinster2012:textual}{divent} for other orders.
#' The result is identical to the total entropy calculated by
#' [ent_phylo].
#' It is much faster but no bias correction is available.
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
#' @returns A tibble with the site names, the estimators used and the estimated entropy.
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' # entropy of each community
#' ent_allen(paracou_6_abd, tree = paracou_6_taxo)
#'
#' # gamma entropy
#' ent_allen(paracou_6_abd, tree = paracou_6_taxo, gamma = TRUE)
#'
#' @name ent_allen
NULL


#' @rdname ent_allen
#'
#' @export
ent_allen <- function(x, tree, ...) {
  UseMethod("ent_allen")
}


#' @rdname ent_allen
#'
#' @export
ent_allen.numeric <- function(
    x,
    tree,
    q = 1,
    normalize = TRUE,
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
  # Normalize x to have probabilities
  the_prob <- x / sum(x)
  # Get unnormalized probabilities p(b)
  ltips <- lapply(tree$phylo$edge[, 2], FUN = function(node) tips(tree, node))
  branches <- unlist(
    lapply(ltips, FUN = function(tips.vector) sum(the_prob[tips.vector]))
  )
  # Calculate Tbar but do not normalize l(b)
  T_bar <- sum(lengths * branches)
  # Eliminate unobserved species (or 0log0 will retrun NaN)
  lengths <- lengths[branches != 0]
  branches <- branches[branches != 0]

  # Sum the lengths of branches with abundance > 0
  the_entropy <- -sum(lengths *  branches^q * ln_q(branches, q)) /
    ifelse(normalize, T_bar, 1)

  # Return
  if (as_numeric) {
    return(the_entropy)
  } else {
    return(
      tibble::tibble_row(
        estimator = "naive",
        order = q,
        entropy = the_entropy
      )
    )
  }
}


#' @rdname ent_allen
#'
#' @export
ent_allen.species_distribution <- function(
    x,
    tree,
    q = 1,
    normalize = TRUE,
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
    the_entropy <- ent_allen.numeric(
      abd,
      tree = tree,
      q = q,
      normalize = normalize,
      as_numeric = as_numeric,
      check_arguments = FALSE
    )
    if (!as_numeric) {
      # Add the site column
      the_entropy <- dplyr::bind_cols(
        site = "Metacommunity",
        the_entropy
      )
    }
    return(the_entropy)
   } else {
    # Apply ent_allen.numeric() to each site
    ent_allen_list <- apply(
      # Eliminate site and weight columns
      x[, !colnames(x) %in% non_species_columns],
      # Apply to each row
      MARGIN = 1,
      FUN = ent_allen.numeric,
      # Arguments
      tree = tree,
      q = q,
      normalize = normalize,
      prune = prune,
      as_numeric = FALSE,
      check_arguments = FALSE
    )
    # Make a tibble with site, estimator and entropy
    the_entropy <- tibble::tibble(
      # Restore non-species columns
      x[colnames(x) %in% non_species_columns],
      # Coerce the list returned by apply into a dataframe
      do.call(rbind.data.frame, ent_allen_list)
    )
    if (as_numeric) {
      return(the_entropy$entropy)
    } else {
      return(the_entropy)
    }
  }
}
