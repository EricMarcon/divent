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
"paracou_6_taxo"



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


#' check_divent_args
#'
#' Checks the arguments of a function of the package divent
#'
#' The function compares the arguments passed to its parent function to the type 
#' they should be and performs some extra tests, *e.g.* probabilities must be positive and sum to 1. 
#' It stops if an argument is not correct.
#' 
#' @param abundances An object of class [abundances].
#' @param alpha The risk level, 5% by default.
#' @param as_numeric If `TRUE`, a number or a numeric vector is returned rather than a tibble.
#' @param check_arguments If `TRUE`, the function arguments are verified.
#' Should be set to `FALSE` to save time when the arguments have been checked elsewhere.
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
#' @param richness_estimator An estimator of richness to evaluate the total number of species,
#' see [div_richness]. Used for interpolation and extrapolation.
#' @param sample_coverage The sample coverage of `x` calculated elsewhere. 
#' Used to calculate the gamma diversity of meta-communities, see details. 
#' @param show_progress If TRUE, a progress bar is shown during long computations. 
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
    warning("Function arguments cannot be checked because the SpatDiv package is not attached. Add CheckArguments=FALSE to suppress this warning or run library('SpatDiv')")
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
  # estimator is checked by match.arg()
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
  # tree
  if (!is.na(names(args["tree"]))) {
    tree <- eval(expression(tree), parent.frame())
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
  # unveiling is checked by match.arg()
  # weights
  if (!is.na(names(args["weights"]))) {
    weights <- eval(expression(weights), parent.frame())
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
  
  # All tests passed.
  return (TRUE)
}


#' Error message
#' 
#' Utility for [check_divent_args]]
#' 
#' @noRd
#'
#' @param message The message to print.
#' @param argument The function argument that did not pass the tests.
#' @param parent_function The name of the function the argument was passed to.
#'
#' @return Nothing
error_message <- function(message, argument, parent_function) {
  cat(deparse(substitute(argument)), "cannot be:\n")
  print(utils::head(argument))
  cat(paste(paste("Error in ", parent_function, ":"), message, "\n"))
  stop("Check the function arguments.", call. = FALSE)
}


#' Check values of a vector are integers
#' 
#' Check that the values of a vector are integer, whatever their type.
#' 
#' @noRd
#'
#' @param x A numeric vector.
#'
#' @return `TRUE` if values are integers.
#'
is_integer_values <- function (x) {
  x_int <- round(x)
  # Return TRUE if no value in x has been modified by rounding
  return(!(any(abs(x_int - x) > sum(x) * .Machine$double.eps)))
}
