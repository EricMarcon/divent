#  Package description ----
#' divent
#'
#' Measures of Diversity and Entropy
#' 
#' This package is a reboot of the **entropart** package \insertCite{Marcon2014c}{divent}.
#'
#' @name divent
#' @docType package
#' @importFrom Rdpack reprompt
#' 
#' @references
#' \insertAllCited{}
NULL

#  Initialization ----
# Columns to ignore when computing species distributions 
non_species_columns <- c(
  # Site characteristics
  "site",
  "weight",
  # Phylodiversity cuts
  "cut",
  "interval",
  # Outputs
  "q",
  "entropy",
  "diversity",
  "abundance",
  # Free comments
  "comments"
)



#  Data ----
#' Paracou plot 6
#'
#' A community assembly.
#' It contains number of trees per species of the plot #6 of Paracou.
#' The plot covers 6.25 ha of tropical rainforest, divided into 4 equally-sized subplots.
#' Each line of the tibble is a subplot. 
#' The "site" column contains the subplot number, "weight" contains its area and all others columns contain a species.
#' Data are the number of trees above 10 cm diameter at breast height (DBH).
#' 
#' This dataset is from Paracou field station, French Guiana, managed by [Cirad](https://www.cirad.fr).
#'
#' @format An object of class [abundances], which is also a [tibble::tibble].
#' @source Permanent data census of Paracou: <https://paracou.cirad.fr/>
#' @seealso [paracou_6_taxo], [paracou_6_fundist]
#' 
"paracou_6_abd"

#' Taxonomy of Paracou plot 6 species
#'
#' The taxonomy of species of the dataset [paracou_6_abd].
#' Distances in the tree are 1 (different species of the same genus),
#' 2 (same family) or 3 (different families).
#' 
#' This dataset is from Paracou field station, French Guiana, managed by [Cirad](https://www.cirad.fr).
#'
#' @format An object of class [ape::phylo], which is a phylogenetic tree.
#' @source Permanent data census of Paracou: <https://paracou.cirad.fr/>
#' @seealso [paracou_6_abd], [paracou_6_fundist]
#' 
"paracou_6_taxo"


#' Functional distances between Paracou plot 6 species
#'
#' A functional distance matrix of species of the dataset [paracou_6_abd].
#' Distances were computed from a trait dataset including specific leaf area, 
#' wood density, seed mass and 95th percentile of height. Gower's metric
#' \insertCite{Gower1971}{divent} was used to obtain a distance matrix.
#' 
#' This dataset is from Paracou field station, French Guiana, managed by [Cirad](https://www.cirad.fr).
#'
#' @format A matrix.
#' @source Permanent data census of Paracou: <https://paracou.cirad.fr/>
#' @seealso [paracou_6_abd], [paracou_6_taxo]
#' 
"paracou_6_fundist"


#' Mock data
#'
#' A simple dataset to test diversity functions.
#' It contains 3 species with their abundances, their distance matrix and
#' their phylogenetic tree.
#'
#' @format `mock_3sp_abd` is a vector, `mock_3sp_dist` a matrix and 
#' `mock_3sp_tree` an object of class [ape::phylo].
#'
#' @name mock_3sp
#' @examples 
#' mock_3sp_abd
#' mock_3sp_dist
#' plot(mock_3sp_tree)
#' axis(2)
#' 
"mock_3sp_abd"

#' @rdname  mock_3sp
"mock_3sp_dist"

#' @rdname  mock_3sp
"mock_3sp_tree"



#  Utilities ----

# Names of variables inside functions:
# abd: a numeric vector of abundances
# prob: a numeric vector of probabilities
# prob_unv : unveiled probabilities
# abundances / probabilities: an object of class abundances / probabilities
# s_0, s_1,  ...: species observed 0, 1, ... times
# s_obs: number of observed species
# s_est: estimated number of species (asymptotic)
# sample_size: number of observed individuals
# sample_coverage: coverage of order 1
# coverage_deficit_2: coverage deficit of order 2 of the sample
# species_names: char vector with the species names of the community


#' Chao's A
#' 
#' Helper for Chao's estimators.
#' 
#' A's formula \insertCite{@Chao2015@, eq. 6b}{divent}) depends on the presence
#' of singletons and doubletons.
#' 
#' @param abd A vector of positive integers.
#'
#' @return The value of A.
#' @noRd
#'
chao_A <- function(abd) {
  
  # Calculate abundance distribution
  abd_distribution <- tapply(abd, INDEX = abd, FUN = length)
  s_1 <- as.numeric(abd_distribution["1"])
  s_2 <- as.numeric(abd_distribution["2"])
  sample_size <- sum(abd)
  
  # Calculate A
  if (is.na(s_1)) {
    A <- 0
  } else {
    # Use Chao1 estimator to evaluate the number of unobserved species
    if (is.na(s_2)) {
      s_0 <- (sample_size - 1) * s_1 * (s_1 - 1) / 2 / sample_size
    } else {
      s_0 <- (sample_size - 1) * s_1^2 / 2 / sample_size / s_2
    }
    A <- 1 - sample_size * s_0 / (sample_size * s_0 + s_1)
  }
  
  return(A)
}


#' check_divent_args
#'
#' Checks the arguments of a function of the package divent
#'
#' The function compares the arguments passed to its parent function to the type 
#' they should be and performs some extra tests, *e.g.* probabilities must be positive and sum to 1. 
#' It stops if an argument is not correct.
#' 
#' The function is always called without arguments.
#' Its arguments exist only for documentation.
#' 
#' @param abundances An object of class [abundances].
#' @param alpha The risk level, 5% by default.
#' @param as_numeric If `TRUE`, a number or a numeric vector is returned rather than a tibble.
#' @param check_arguments If `TRUE`, the function arguments are verified.
#' Should be set to `FALSE` to save time when the arguments have been checked elsewhere.
#' @param coverage_estimator An estimator of sample coverage used by [coverage].
#' @param distances A distance matrix or an object of class [stats::dist]
#' @param estimator An estimator of asymptotic entropy, diversity or richness.
#' @param gamma If `TRUE`, \eqn{\gamma} diversity, i.e. diversity of the metacommunity, is computed.
#' @param jack_alpha The risk level, 5% by default, used to optimize the jackknife order.
#' @param jack_max The highest jackknife order allowed. Default is 10.
#' @param level The level of interpolation or extrapolation. 
#' It may be a sample size (an integer) or a sample coverage 
#' (a number between 0 and 1).
#' If not `NULL`, the asymptotic `estimator` is ignored.
#' @param n_simulations The number of simulations used to estimate the confidence envelope.
#' @param normalize If `TRUE`, phylogenetic is normalized: the height of the tree is set to 1.
#' @param probability_estimator A string containing one of the possible estimators
#' of the probability distribution (see [probabilities]). 
#' Used only for extrapolation.
#' @param q The order of diversity.
#' @param rate The decay rate of the exponential similarity.
#' @param richness_estimator An estimator of richness to evaluate the total number of species,
#' see [div_richness]. Used for interpolation and extrapolation.
#' @param sample_coverage The sample coverage of `x` calculated elsewhere. 
#' Used to calculate the gamma diversity of meta-communities, see details. 
#' @param show_progress If TRUE, a progress bar is shown during long computations. 
#' @param similarities A similarity matrix, that can be obtained by [fun_similarity].
#' Its default value is the identity matrix.
#' @param species_distribution An object of class [species_distribution].
#' @param tree An ultrametric, phylogenetic tree.
#' May be an object of class [phylo_divent], [ape::phylo], [ade4::phylog] or [stats::hclust]. 
#' @param unveiling A string containing one of the possible unveiling methods 
#' to estimate the probabilities of the unobserved species (see [probabilities]).
#' Used only for extrapolation.
#' @param weights The weights of the sites of the species distributions.
#'
#' @return Returns `TRUE` or stops if a problem is detected.
#' 
#' @export
#'
#' @keywords internal
#' 
check_divent_args <- function(
    abundances = NULL,
    as_numeric = NULL,
    check_arguments = NULL,
    estimator = NULL,
    jack_alpha = NULL,
    jack_max = NULL,
    level = NULL,
    n_simulations = NULL,
    probability_estimator = NULL,
    q = NULL,
    richness_estimator = NULL,
    sample_coverage = NULL,
    show_progress = NULL,
    unveiling = NULL) {

  # Verify that the package is attached
  if (!"divent" %in% .packages()) {
    warning("Function arguments cannot be checked because the package divent is not attached. Add CheckArguments=FALSE to suppress this warning or run library('SpatDiv')")
    return (TRUE)
  }
  # Get the list of arguments of the parent function
  parent_function <- sys.call(-1)[[1]]
  # If apply() or similar was used, the function name is not in parent_function: sys.call(-1)[[1]] returns "FUN"
  if (parent_function == "FUN") {
    warning("Function arguments cannot be checked, probably because you used apply(). Add CheckArguments=FALSE to suppress this warning.")
    return (TRUE)
  }

  # Find the arguments. match.fun does not work with package::function
  # as.character creates a vector. The name of the function is the last item
  parent_function_split <- as.character(parent_function)
  parent_function_name <- parent_function_split[length(parent_function_split)]
  args <- formals(match.fun(parent_function_name))
  
  # abundances
  if (!is.na(names(args["abundances"]))) {
    abundances <- eval(expression(abundances), parent.frame())
    if (!is_abundances(abundances)) {
      error_message(
        "abundances must be an object of class 'abundances'", 
        abundances, 
        parent_function
      )
    }
  }
  # alpha
  if (!is.na(names(args["alpha"]))) {
    alpha <- eval(expression(alpha), parent.frame())
    if (!is.numeric(alpha) | length(alpha)!=1) {
      error_message(
        "alpha must be a number.",
        alpha,
        parent_function
      )
    }
    if (any(alpha < 0) | any(alpha > 1)) {
      error_message(
        "alpha must be between 0 and 1",
        alpha,
        parent_function
      )
    }
  }
  # as_numeric
  if (!is.na(names(args["as_numeric"]))) {
    as_numeric <- eval(expression(as_numeric), parent.frame())
    if (!is.logical(as_numeric) | length(as_numeric) != 1) {
      error_message(
        "as_numeric must be TRUE or FALSE", 
        as_numeric, 
        parent_function
      )
    }
  }
  # check_arguments
  if (!is.na(names(args["check_arguments"]))) {
    check_arguments <- eval(expression(check_arguments), parent.frame())
    if (!is.logical(check_arguments) | length(check_arguments) != 1) {
      error_message(
        "check_arguments must be TRUE or FALSE", 
        check_arguments, 
        parent_function
      )
    }
  }
  # coverage_estimator is checked by match.arg()
  # estimator is checked by match.arg()
  # distances
  if (!is.na(names(args["distances"]))) {
    distances <- eval(expression(distances), parent.frame())
    if (!is.null(distances) && !inherits(distances, "dist")) {
      if (!is.numeric(distances)) {
        error_message(
          "distances must be a numeric",
          distances,
          parent_function
        )
      }
      if (!is.matrix(distances)) {
        error_message(
          "distances must be a matrix",
          distances,
          parent_function
        )
      }
      if (nrow(distances) != ncol(distances)) {
        error_message(
          "distances must be a square matrix",
          distances,
          parent_function
        )
      }
      if (!isSymmetric(distances)) {
        error_message(
          "distances must be a symmetric matrix",
          distances,
          parent_function
        )
      }
      if (any(distances < 0)) {
        error_message(
          "distances must be positive",
          distances,
          parent_function
        )
      }
      if (any(diag(distances != 0))) {
        error_message(
          "distances must be zero between a species and itself",
          distances,
          parent_function
        )
      }
    }
  }
  # gamma
  if (!is.na(names(args["gamma"]))) {
    gamma <- eval(expression(gamma), parent.frame())
    if (!is.logical(gamma) | length(gamma) != 1) {
      error_message(
        "gamma must be TRUE or FALSE", 
        gamma, 
        parent_function
      )
    }
  }
  # jack_alpha
  if (!is.na(names(args["jack_alpha"]))) {
    jack_alpha <- eval(expression(jack_alpha), parent.frame())
    if (!is.numeric(jack_alpha) | length(jack_alpha)!=1) {
      error_message(
        "jack_alpha must be a number.",
        jack_alpha,
        parent_function
      )
    }
    if (any(jack_alpha <= 0) | any(jack_alpha >= 1)) {
      error_message(
        "jack_alpha must be between 0 and 1",
        jack_alpha,
        parent_function
      )
    }
  }
  # jack_max
  if (!is.na(names(args["jack_max"]))) {
    jack_max <- eval(expression(jack_max), parent.frame())
    if (!is.numeric(jack_max) | length(jack_max)!=1) {
      error_message(
        "jack_max must be a number.",
        jack_max,
        parent_function
      )
    }
    if (any(jack_max < 1) | any(jack_max > 10)) {
      error_message(
        "jack_max must be between 1 and 10",
        jack_max,
        parent_function
      )
    }
  }
  # level
  if (!is.na(names(args["level"]))) {
    level <- eval(expression(level), parent.frame())
    if (!is.null(level)) {
      if (!is.numeric(level) | length(level)!=1) {
        error_message(
          "level must be a number.",
          level,
          parent_function
        )
      }
      if (any(level <=0)) {
        error_message(
          "level must be positive.",
          level,
          parent_function
        )
      }
    }
   }
  # n_simulations
  if (!is.na(names(args["n_simulations"]))) {
    n_simulations <- eval(expression(n_simulations), parent.frame())
    if (!is.numeric(n_simulations) | length(n_simulations)!=1) {
      error_message(
        "n_simulations must be a number.",
        n_simulations,
        parent_function
      )
    }
    if (any(n_simulations !=0 & n_simulations < 2)) {
      error_message(
        "n_simulations must be 0 or at least 2",
        n_simulations,
        parent_function
      )
    }
  }
  # normalize
  if (!is.na(names(args["normalize"]))) {
    normalize <- eval(expression(normalize), parent.frame())
    if (!is.logical(normalize) | length(normalize) != 1) {
      error_message(
        "normalize must be TRUE or FALSE", 
        normalize, 
        parent_function
      )
    }
  }
  # probability_estimator is checked by match.arg()
  # richness_estimator is checked by match.arg()
  # q
  if (!is.na(names(args["q"]))) {
    q <- eval(expression(q), parent.frame())
    if (!is.numeric(q) | length(q) != 1) {
      error_message(
        "q must be a number.", 
        q, 
        parent_function
      )
    }
  }
  # rate
  if (!is.na(names(args["rate"]))) {
    rate <- eval(expression(rate), parent.frame())
    if (!is.null(rate)) {
      if (!is.numeric(rate) | length(rate)!=1) {
        error_message(
          "rate must be a number.",
          rate,
          parent_function
        )
      }
      if (any(rate <=0)) {
        error_message(
          "rate must be positive.",
          rate,
          parent_function
        )
      }
    }
  }
  # sample_coverage
  if (!is.na(names(args["sample_coverage"]))) {
    sample_coverage <- eval(expression(sample_coverage), parent.frame())
    if (!is.null(sample_coverage)) {
      if (!is.numeric(sample_coverage) | length(sample_coverage)!=1) {
        error_message(
          "sample_coverage must be a number.",
          sample_coverage,
          parent_function
        )
      }
      if (any(sample_coverage <= 0) | any(sample_coverage >= 1)) {
        error_message(
          "sample_coverage must be between 0 and 1",
          sample_coverage,
          parent_function
        )
      }
    }
  }
  # show_progress
  if (!is.na(names(args["show_progress"]))) {
    show_progress <- eval(expression(show_progress), parent.frame())
    if (!is.logical(show_progress) | length(show_progress) != 1) {
      error_message(
        "show_progress must be TRUE or FALSE", 
        show_progress, 
        parent_function
      )
    }
  }
  # similarities
  if (!is.na(names(args["similarities"]))) {
    similarities <- eval(expression(similarities), parent.frame())
    if (!is.numeric(similarities)) {
      error_message(
        "similarities must be a numeric",
        similarities,
        parent_function
      )
    }
    if (!is.matrix(similarities)) {
      error_message(
        "similarities must be a matrix",
        similarities,
        parent_function
      )
    }
    if (nrow(similarities) != ncol(similarities)) {
      error_message(
        "similarities must be a square matrix",
        similarities,
        parent_function
      )
    }
    if (any(similarities = 0 | similarities > 1)) {
      error_message(
        "similarities must be between 0 and 1",
        similarities,
        parent_function
      )
    }
    if (any(diag(similarities != 1))) {
      error_message(
        "similarities must be 1 between a species and itself",
        similarities,
        parent_function
      )
    }
  }
  # species_distribution
  if (!is.na(names(args["species_distribution"]))) {
    species_distribution <- eval(expression(species_distribution), parent.frame())
    if (!is_species_distribution(species_distribution)) {
      error_message(
        "species_distribution must be an object of class 'species_distribution'", 
        species_distribution, 
        parent_function
      )
    }
  }
  # tree
  if (!is.na(names(args["tree"]))) {
    tree <- eval(expression(tree), parent.frame())
    if (!is.null(tree)) {
      if (
          !inherits(tree, "phylo_divent") &
          !inherits(tree, "phylo") &
          !inherits(tree, "phylog") &
          !inherits(tree, "hclust")) {
        error_message(
          "tree must be an object of class 'phylo_divent', 'phylo', 'phylog' or 'hclust'", 
          tree, 
          parent_function
        )
      }
    }
  }
  # unveiling is checked by match.arg()
  # weights
  if (!is.na(names(args["weights"]))) {
    weights <- eval(expression(weights), parent.frame())
    if (!is.null(weights)) {
      if (!is.numeric(weights)) {
        error_message(
          "weights must be a numeric vector",
          weights,
          parent_function
        )
      }
      if (any(weights < 0)) {
        error_message(
          "weights must be positive",
          weights,
          parent_function
        )
      }
    }
  }
  
  # All tests passed.
  return (TRUE)
}


#' Check Similarity or Distance Matrix
#' 
#' Verify that a similarity or distance matrix fits a species distribution and 
#' filter it so that its elements are the same, in the same order.
#' 
#' If species names are missing, just check the dimensions.
#'
#' @param sim_dist_matrix A similarity or distance matrix, or a [stats::dist].
#' @param species_distribution A species distribution, or a named vector.
#'
#' @return A similarity matrix that corresponds to the species distribution.
#' @noRd
#'
checked_matrix <- function(
    sim_dist_matrix,
    species_distribution) {
  
  if (inherits(sim_dist_matrix, "dist")) {
    # dist objects are supported but the remainder assumes a matrix
    sim_dist_matrix <- as.matrix(sim_dist_matrix)
  }
  
  # No names needed
  if (
    # No names in the matrix
    is.null(colnames(sim_dist_matrix)) | 
    # No names in the distribution that may be a tibble or a vector
    (is.null(colnames(species_distribution)) & is.null(names(species_distribution)))
  ) {
    # The matrix may not be named
    if (ncol(sim_dist_matrix) == length(species_distribution)) {
      # Do not change the matrix
      return(sim_dist_matrix)
    } else {
      stop("If the similarity matrix is not named, then its size must fit the number of species.")
    }
  }
  
  # Get species names
  if (is_species_distribution(species_distribution)) {
    is_species_column <- !colnames(species_distribution) %in% non_species_columns
    species_names <- colnames(species_distribution)[is_species_column]
  } else if (is.vector(species_distribution)) {
    species_names <- names(species_distribution)
  }
  
  # Stop if some species are not in the matrix
  if (length(species_names) == 0) stop("There are no species in the distribution")
  if (length(setdiff(species_names, colnames(sim_dist_matrix))) != 0) {
    stop("Some species are missing in the similarity matrix.")    
  } 
  
  # Filter and reorder the similarity matrix
  return(sim_dist_matrix[species_names, species_names])
}


#' Gamma entropy of a metacommunity
#' 
#' Build the metacommunity and check that abundances are integers.
#' 
#' If they are not (due to weights of communities) then use a fallback estimator:
#' "ChaoShen" requires the sample coverage of the assemblage of sites.
#' "Grassberger" accepts non integer abundances.
#' "Marcon" combines both.
#' 
#' [ent_tsallis.numeric] contains the only implementation of this estimation.
#' i.e., [ent_shannon] can't be used but [ent_tsallis.numeric] 
#' with `q=1` will work fine.
#' 
#' @param species_distribution An object of class [species_distribution].
#'
#' @return A tibble with the estimator used and the estimated entropy.
#' @noRd
#' 
ent_gamma_hill <- function(
    species_distribution,
    q,
    estimator,
    level,
    probability_estimator,
    unveiling,
    richness_estimator,
    jack_alpha,
    jack_max,
    coverage_estimator,
    as_numeric) {
  
  # Build the metacommunity
  abd <- metacommunity.abundances(
    species_distribution, 
    as_numeric = TRUE, 
    check_arguments = FALSE
  )
  if (is_integer_values(abd)) {
    # Sample coverage is useless
    sample_coverage <- NULL
  } else {
    # Non-integer values in the metacommunity. 
    # Calculate the sample coverage and change the estimator.
    sample_coverage <- coverage.numeric(
      colSums(
        species_distribution[
          , !colnames(species_distribution) %in% non_species_columns
        ]
      ),
      estimator = coverage_estimator,
      as_numeric = TRUE,
      check_arguments = FALSE
    )
    if (!estimator %in% c("Marcon", "ChaoShen")) {
      estimator <- "Marcon"
    }
  }
  
  # Compute the entropy. Call the appropriate function for its estimators.
  # Richness estimators are specific
  if (q==0 & estimator %in% c("jackknife", "iChao1", "Chao1", "rarefy", "naive")) {
    the_diversity <- div_richness.numeric(
      abd, 
      estimator = estimator,
      jack_alpha  = jack_alpha,
      jack_max = jack_max,
      level = level, 
      probability_estimator = probability_estimator,
      unveiling = unveiling,
      coverage_estimator = coverage_estimator,
      as_numeric = FALSE,
      check_arguments = FALSE
    )
    # Calculate entropy
    the_entropy <- dplyr::mutate(
      the_diversity, 
      entropy = .data$diversity - 1, 
      .keep = "unused"
    )
    # entropy must be a number if as_numeric = TRUE
    if (as_numeric) {
      the_entropy <- the_entropy$entropy
    }
  } else if (q==1 & is.null(sample_coverage)) {
    # Non-integer values in the metacommunity are supported only by ent_tsallis
    the_entropy <- ent_shannon.numeric(
      abd,
      estimator = estimator,
      level = level, 
      probability_estimator = probability_estimator,
      unveiling = unveiling,
      richness_estimator = richness_estimator,
      jack_alpha  = jack_alpha,
      jack_max = jack_max,
      coverage_estimator = coverage_estimator,
      as_numeric = as_numeric,
      check_arguments = FALSE
    )
  } else if (q==2 & is.null(sample_coverage)) {
    # Non-integer values in the metacommunity are supported only by ent_tsallis
    the_entropy <- ent_simpson.numeric(
      abd,
      estimator = estimator,
      level = level, 
      probability_estimator = probability_estimator,
      unveiling = unveiling,
      richness_estimator = richness_estimator,
      jack_alpha  = jack_alpha,
      jack_max = jack_max,
      coverage_estimator = coverage_estimator,
      as_numeric = as_numeric,
      check_arguments = FALSE
    )
  } else {
    the_entropy <- ent_tsallis.numeric(
      abd,
      q = q,
      estimator = estimator,
      level = level, 
      probability_estimator = probability_estimator,
      unveiling = unveiling,
      richness_estimator = richness_estimator,
      jack_alpha  = jack_alpha,
      jack_max = jack_max,
      sample_coverage = sample_coverage,
      as_numeric = as_numeric,
      check_arguments = FALSE
    )
  }
  # Add the site column
  if (!as_numeric) {
    the_entropy <- dplyr::bind_cols(
      site = "Metacommunity",
      the_entropy
    )
  }
  return(the_entropy)
}


#' Similarity-Based Gamma entropy of a metacommunity
#' 
#' Build the metacommunity and check that abundances are integers.
#' 
#' See [ent_gamma_hill] for details.
#' @param species_distribution An object of class [species_distribution].
#' 
#' @return A tibble with the estimator used and the estimated entropy.
#' @noRd
#' 
ent_gamma_similarity <- function(
    species_distribution,
    similarities,
    q,
    estimator,
    probability_estimator,
    unveiling,
    jack_alpha,
    jack_max,
    coverage_estimator,
    as_numeric) {
  
  # Build the metacommunity
  abd <- metacommunity.abundances(
    species_distribution, 
    as_numeric = TRUE, 
    check_arguments = FALSE
  )
  if (is_integer_values(abd)) {
    # Sample coverage is useless
    sample_coverage <- NULL
  } else {
    # Non-integer values in the metacommunity. 
    # Calculate the sample coverage and change the estimator.
    sample_coverage <- coverage.numeric(
      colSums(
        species_distribution[
          , !colnames(species_distribution) %in% non_species_columns
        ]
      ),
      estimator = coverage_estimator,
      as_numeric = TRUE,
      check_arguments = FALSE
    )
    if (!estimator %in% c("Marcon", "ChaoShen")) {
      estimator <- "Marcon"
    }
  }
  
  # Compute the entropy.
  the_entropy <- ent_similarity.numeric(
    abd,
    similarities = similarities,
    q = q,
    estimator = estimator,
    probability_estimator = probability_estimator,
    unveiling = unveiling,
    jack_alpha  = jack_alpha,
    jack_max = jack_max,
    sample_coverage = sample_coverage,
    as_numeric = as_numeric,
    check_arguments = FALSE
  )
  
  # Return
  if (as_numeric) {
    return(the_entropy)
  } else {
    return(
      dplyr::bind_cols(
        site = "Metacommunity",
        the_entropy
      )
    )
  }
}


#' Error message
#' 
#' Utility for [check_divent_args]
#' 
#' @param message The message to print.
#' @param argument The function argument that did not pass the tests.
#' @param parent_function The name of the function the argument was passed to.
#'
#' @return Nothing
#' @noRd
#' 
error_message <- function(message, argument, parent_function) {
  cat(deparse(substitute(argument)), "cannot be:\n")
  print(utils::head(argument))
  cat(paste(paste("Error in ", parent_function, ":"), message, "\n"))
  stop("Check the function arguments.", call. = FALSE)
}


#' Check Integers
#' 
#' Check that the values of a vector are integer, whatever their type.
#' 
#' @param x A numeric vector.
#'
#' @return `TRUE` if values are integers.
#' @noRd
#'
is_integer_values <- function(x) {
  x_int <- round(x)
  # Return TRUE if no value in x has been modified by rounding
  return(!any(abs(x_int - x) > sum(x) * .Machine$double.eps))
}


#' Phylogenetic abundances
#' 
#' Calculate abundances of species and their ancestors along a phylogenetic tree.
#' 
#' @return A list of matrices. 
#' Each matrix contains the abundances of species (in lines) of each community
#' (in columns) in an interval of the tree.
#' @noRd
#'
phylo_abd <- function(
    abundances,
    tree) {

  # Calculate abundances along the tree, that are a list of matrices
  sapply(
    # Each phylogenetic group yields an item of the list
    colnames(tree$phylo_groups), 
    function(group) {
      # Create a matrix with the abundances of groups in each community
      apply(
        abundances[, !colnames(abundances) %in% non_species_columns], 
        # Each community yields a column of the matrix
        MARGIN = 1, 
        FUN = function(abd) {
          tapply(
            as.numeric(abd), 
            # Each group yields a row of the matrix
            INDEX = tree$phylo_groups[names(abd), group], 
            FUN = sum
          )
        }
      )
    }, 
    simplify = FALSE
  )
}
  

#' Phylogenetic entropies
#' 
#' Calculate entropies of a list of phylogenetic abundances (obtained by 
#' [phylo_abd]).
#' Each item of the list corresponds to a phylogenetic group, i.e. an interval
#' of the tree (where the species do not change).
#' 
#' @param phylo_abd A list of matrices of abundance (caution: lines are species,
#' columns are communities).
#'
#' @return A vector
#'  
#' @noRd
#'
phylo_entropy.phylo_abd <- function(
    # Allow lapply along q
    q,
    phylo_abd,
    tree,
    normalize,
    # Other arguments for ent_tsallis.numeric
    estimator,
    level, 
    probability_estimator,
    unveiling,
    richness_estimator,
    jack_alpha, 
    jack_max,
    coverage_estimator) {
 
  # Calculate entropy of each community in each group.
  # simplify2array() makes a matrix with the list of vectors.
  phylo_entropies <- simplify2array(
    lapply(
      # Calculate entropy in each item of the list, i.e. group.
      # Obtain a list.
      phylo_abd,
      FUN = function(group) {
        apply(
          group,
          # Calculate entropy of each column of the matrix, i.e. community.
          MARGIN = 2,
          FUN = ent_tsallis.numeric,
          # Arguments
          q = q,
          estimator = estimator,
          level = level, 
          probability_estimator = probability_estimator,
          unveiling = unveiling,
          richness_estimator = richness_estimator,
          jack_alpha  = jack_alpha, 
          jack_max = jack_max,
          coverage_estimator = coverage_estimator,
          # Obtain a vector.
          as_numeric = TRUE,
          check_arguments = FALSE
        )
      }
    )
  )
  # Should be a matrix, but simplify2array() makes a vector instead of a 1-col 
  # matrix. Force a matrix.
  if(is.vector(phylo_entropies)) {
    phylo_entropies <- t(phylo_entropies)
  }
  
  # Calculate the weighted mean of entropy and normalize
  the_entropy <- as.numeric(tree$intervals %*% t(phylo_entropies))
  if (normalize) the_entropy <- the_entropy / sum(tree$intervals)
  
  return(the_entropy)
}
