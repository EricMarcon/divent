# Package description ----
#' divent
#'
#' Measures of Diversity and Entropy
#'
#' This package is a reboot of the **entropart** package \insertCite{Marcon2014c}{divent}.
#'
#' @importFrom Rdpack reprompt
#'
#' @references
#' \insertAllCited{}
"_PACKAGE"


# C code ----
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppParallel RcppParallelLibs
#' @useDynLib divent, .registration = TRUE


# Functions to reexport ----
#' @export
ggplot2::autoplot
base::as.numeric
base::as.matrix


#' ggplot method to plot wmppp objects
#'
#' This method is from the dbmss package. See [dbmss::autoplot.wmppp].
#'
#' @param object an object to be plotted.
#' @param ... extra arguments, currently unused.
#' @param show.window if `TRUE`, the borders of the window containing the points are shown on the point map.
#' @param MaxPointTypes the maximum number of different point types to show.
#' If the point set contains more of them, the less frequent ones are gathered as "Other".
#' This number must be limited for readability and not to exceed the number of colors offered by the palette.
#' @param Other the name of the point types gathered as "Other"
#' @param main the title of the plot.
#' @param xlab the X-axis label.
#' @param ylab the Y-axis label.
#' @param LegendLabels a vector of characters.
#' The first two items describe the observed and null-hypothesis curves, the third and last item the confidence interval.
#' To be used only in plots with two curves (typically observed and expected values).
#' The default is `NULL` to display the full description of functions.
#' @param labelSize the guide of the point size legend in point maps, i.e. what the `PointSize` mark represents.
#' @param labelColor the guide of the point color legend in point maps, i.e. what the `PointType` mark represents.
#' @param palette The color palette used to display point types in maps. See [ggplot2::scale_colour_brewer].
#' @param windowColor the color used to draw the limits of the windows in point maps.
#' @param windowFill the color used to fill the windows in point maps.
#' @param alpha the opacity of the confidence envelope (in function values) or the points (in maps), between 0 and 1.
#'
#' @returns A [ggplot2::ggplot].
#' @export
#'
#' @examples
#' autoplot(paracou_6_wmppp)
#'
autoplot.wmppp <- function(
    object,
    ...,
    show.window = TRUE,
    MaxPointTypes = 6,
    Other = "Other",
    main = NULL,
    xlab = NULL,
    ylab = NULL,
    LegendLabels = NULL,
    labelSize = "Weight",
    labelColor = "Type",
    palette="Set1",
    windowColor = "black",
    windowFill = "transparent",
    alpha = 1) {

  return(
    dbmss:::autoplot.wmppp(
      object = object,
      ...,
      show.window = show.window,
      MaxPointTypes = 6,
      Other = Other,
      main = main,
      xlab = xlab,
      ylab = ylab,
      LegendLabels = LegendLabels,
      labelSize = labelSize,
      labelColor = labelColor,
      palette = palette,
      windowColor = windowColor,
      windowFill = windowFill,
      alpha = alpha
    )
  )
}


# Initialization ----

utils::globalVariables("non_species_columns")


# Data ----

#' Non-Species columns
#'
#' The column names that are not those of species in a [species_distribution].
#'
#' Objects of classes [abundances] and [probabilities], that are also of class
#' [species_distribution], have columns named after the species they contain.
#' Some columns are reserved to describe the plots and their diversity.
#'
#'
#' @format A character vector.
#' @name non_species_columns
"non_species_columns"


#' Paracou plot 6
#'
#' A community assembly.
#' It contains number of trees per species of the plot #6 of Paracou.
#' The plot covers 6.25 ha of tropical rainforest, divided into 4 equally-sized subplots.
#'
#' In `paracou_6_abd` (a tibble), the "site" column contains the subplot number, "weight" contains its area and all others columns contain a species.
#' Data are the number of trees above 10 cm diameter at breast height (DBH).
#'
#' In `paracou_6_wmppp` (a point pattern), the point type is tree species and the point weight is their basal area, in square centimeters.
#'
#' This dataset is from Paracou field station, French Guiana, managed by [Cirad](https://www.cirad.fr).
#'
#' @format `paracou_6_abd` is an object of class [abundances], which is also a [tibble::tibble].
#' Each line of the tibble is a subplot.
#' `paracou_6_wmppp` is a [dbmss::wmppp] object, i.e. a weighted, marked planar point pattern.
#' @source Permanent data census of Paracou: <https://paracou.cirad.fr/>
#' @seealso [paracou_6_taxo], [paracou_6_fundist]
#' @name paracou_6
"paracou_6_abd"

#' @rdname paracou_6
"paracou_6_wmppp"


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
#' @references
#' \insertAllCited{}
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



# Utilities ----

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
#' The function is always called without arguments.
#' Its arguments exist only for documentation.
#'
#' @param abd a numeric vector containing abundances.
#' @param abundances an object of class [abundances].
#' @param alpha the risk level, 5% by default.
#' @param as_numeric if `TRUE`, a number or a numeric vector is returned rather than a tibble.
#' @param bootstrap The method used to obtain the probabilities to generate
#' bootstrapped communities from observed abundances.
#' If "Marcon2012", the probabilities are simply the abundances divided by the total
#' number of individuals \insertCite{Marcon2012a}{divent}.
#' If "Chao2013" or "Chao2015" (by default), a more sophisticated approach is used
#' (see [as_probabilities]) following \insertCite{Chao2013;textual}{divent} or
#' \insertCite{Chao2015;textual}{divent}.
#' @param check_arguments if `TRUE`, the function arguments are verified.
#' Should be set to `FALSE` to save time when the arguments have been checked elsewhere.
#' @param correction the edge-effect correction to apply when estimating
#' the number of neighbors.
#' @param coverage_estimator an estimator of sample coverage used by [coverage].
#' @param distances a distance matrix or an object of class [stats::dist].
#' @param distribution The distribution of species abundances.
#' May be "lnorm" (log-normal), "lseries" (log-series), "geom" (geometric) or
#' "bstick" (broken stick).
#' @param estimator an estimator of asymptotic entropy, diversity or richness.
#' @param fisher_alpha Fisher's \eqn{\alpha} in the log-series distribution.
#' @param gamma if `TRUE`, \eqn{\gamma} diversity, i.e. diversity of the metacommunity, is computed.
#' @param global if `TRUE`, a global envelope sensu \insertCite{Duranton2005}{divent} is calculated.
#' @param jack_alpha the risk level, 5% by default, used to optimize the jackknife order.
#' @param jack_max the highest jackknife order allowed. Default is 10.
#' @param k the order of Hurlbert's diversity.
#' @param level the level of interpolation or extrapolation.
#' It may be a sample size (an integer) or a sample coverage
#' (a number between 0 and 1).
#' If not `NULL`, the asymptotic `estimator` is ignored.
#' @param n the number of observations.
#' @param n_simulations the number of simulations used to estimate the confidence envelope.
#' @param normalize if `TRUE`, phylogenetic is normalized: the height of the tree is set to 1.
#' @param orders The orders of diversity.
#' @param prob a numeric vector containing probabilities.
#' @param prob_geom the proportion of resources taken by successive species
#' of the geometric distribution.
#' @param probability_estimator a string containing one of the possible estimators
#' of the probability distribution (see [probabilities]).
#' Used only for extrapolation.
#' @param q a number: the order of diversity.
#' @param q_threshold the value of `q` above which diversity is computed
#' directly with the naive estimator \eqn{(sum{p_s^q}^{\frac{1}{(1-q)}}},
#' without computing entropy.
#' When `q` is great, the exponential of entropy goes to \eqn{0^{\frac{1}{(1-q)}}},
#' causing rounding errors while the naive estimator of diversity is less and
#' less biased.
#' @param r a vector of distances.
#' @param rate the decay rate of the exponential similarity.
#' @param richness_estimator an estimator of richness to evaluate the total number of species,
#' see [div_richness]. used for interpolation and extrapolation.
#' @param sample_coverage the sample coverage of `x` calculated elsewhere.
#' Used to calculate the gamma diversity of meta-communities, see details.
#' @param show_progress if TRUE, a progress bar is shown during long computations.
#' @param similarities a similarity matrix, that can be obtained by [fun_similarity].
#' Its default value is the identity matrix.
#' @param sd_lnorm the simulated log-normal distribution standard deviation.
#' This is the standard deviation on the log scale.
#' @param size The number of individuals to draw in each community.
#' @param species_number The number of species.
#' @param species_distribution an object of class [species_distribution].
#' @param thomas_scale in Thomas point patterns, the standard deviation of random displacement
#' (along each coordinate axis) of a point from its cluster center.
#' @param thomas_mu in Thomas point patterns, the mean number of points per cluster.
#' The intensity of the Poisson process of cluster centers is calculated as
#' the number of points (`size`) per area divided by `thomas_mu`.
#' @param tree an ultrametric, phylogenetic tree.
#' May be an object of class [phylo_divent], [ape::phylo], [ade4::phylog] or [stats::hclust].
#' @param unveiling a string containing one of the possible unveiling methods
#' to estimate the probabilities of the unobserved species (see [probabilities]).
#' Used only for extrapolation.
#' @param use.names if `TRUE`, the names of the `species_distribution` are kept
#' in the matrix or vector they are converted to.
#' @param weights the weights of the sites of the species distributions.
#' @param w_max the maximum weight in a uniform distribution.
#' @param w_mean the mean weight in an exponential distribution
#' (i.e. the negative of the inverse of the decay rate).
#' @param w_min the minimum weight in a uniform or Weibull distribution.
#' @param weibull_scale the scale parameter in a Weibull distribution.
#' @param weibull_shape the shape parameter in a Weibull distribution.
#' @param X a spatialized community
#' (A [dbmss::wmppp] object with `PointType` values as species names.)
#' @param win the window containing the point pattern.
#' It is an [spatstat.geom::owin] object.
#' Default is a 1x1 square.
#'
#' @returns Returns `TRUE` or stops if a problem is detected.
#'
#' @export
#'
#' @keywords internal
#'
check_divent_args <- function(
    abd = NULL,
    abundances = NULL,
    alpha = NULL,
    as_numeric = NULL,
    bootstrap = NULL,
    check_arguments = NULL,
    correction = NULL,
    coverage_estimator = NULL,
    distances = NULL,
    distribution = NULL,
    estimator = NULL,
    fisher_alpha = NULL,
    gamma = NULL,
    global = NULL,
    jack_alpha = NULL,
    jack_max = NULL,
    k = NULL,
    level = NULL,
    n = NULL,
    n_simulations = NULL,
    normalize = NULL,
    orders = NULL,
    prob = NULL,
    prob_geom = NULL,
    probability_estimator = NULL,
    q = NULL,
    q_threshold = NULL,
    r = NULL,
    rate = NULL,
    richness_estimator = NULL,
    sample_coverage = NULL,
    sd_lnorm = NULL,
    show_progress = NULL,
    similarities = NULL,
    size = NULL,
    species_number = NULL,
    species_distribution = NULL,
    thomas_mu = NULL,
    thomas_scale = NULL,
    tree = NULL,
    use.names = NULL,
    unveiling = NULL,
    weights = NULL,
    w_max = NULL,
    w_mean = NULL,
    w_min = NULL,
    weibull_scale = NULL,
    weibull_shape = NULL,
    X = NULL,
    win = NULL) {

  # Verify that the package is attached
  if (!"divent" %in% .packages()) {
    warning("Function arguments cannot be checked because the package divent is not attached. Add CheckArguments=FALSE to suppress this warning or run library('divent')")
    return(TRUE)
  }
  # Get the list of arguments of the parent function
  parent_function <- sys.call(-1)[[1]]
  # If apply() or similar was used, the function name is not in parent_function: sys.call(-1)[[1]] returns "FUN"
  if (parent_function == "FUN") {
    warning("Function arguments cannot be checked, probably because you used apply(). Add CheckArguments=FALSE to suppress this warning.")
    return(TRUE)
  }

  # Find the arguments. match.fun does not work with package::function
  # as.character creates a vector. The name of the function is the last item
  parent_function_split <- as.character(parent_function)
  parent_function_name <- parent_function_split[length(parent_function_split)]
  args <- formals(match.fun(parent_function_name))

  # abd
  if (!is.na(names(args["abd"]))) {
    abd <- eval(expression(abd), parent.frame())
    if (!is.null(abd)) {
      if (!is.numeric(abd) && !is.vector(abd)) {
        error_message(
          "abd must be a numeric vector",
          abd,
          parent_function
        )
      }
      if (any(abd < 0)) {
        error_message(
          "abundances must not be negative.",
          abd,
          parent_function
        )
      }
      if (!is_integer_values(abd)) {
        error_message(
          "abundances must be integer values.",
          abd,
          parent_function
        )
      }
    }
  }
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
    if (!is.numeric(alpha) | length(alpha) != 1) {
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
  # bootstrap is checked by match.arg()
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
  # correction is checked by match.arg()
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
  # distribution is checked by match.arg()
  # fisher_alpha
  if (!is.na(names(args["fisher_alpha"]))) {
    fisher_alpha <- eval(expression(fisher_alpha), parent.frame())
    if (!is.null(fisher_alpha)) {
      if (!is.numeric(fisher_alpha) | length(fisher_alpha) != 1) {
        error_message(
          "fisher_alpha must be a number.",
          fisher_alpha,
          parent_function
        )
      }
      if (any(fisher_alpha <= 0)) {
        error_message(
          "fisher_alpha must be positive.",
          fisher_alpha,
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
  # global
  if (!is.na(names(args["global"]))) {
    global <- eval(expression(global), parent.frame())
    if (!is.logical(global) | length(global) != 1) {
      error_message(
        "global must be TRUE or FALSE",
        global,
        parent_function
      )
    }
  }
  # jack_alpha
  if (!is.na(names(args["jack_alpha"]))) {
    jack_alpha <- eval(expression(jack_alpha), parent.frame())
    if (!is.numeric(jack_alpha) | length(jack_alpha) != 1) {
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
    if (!is.numeric(jack_max) | length(jack_max) != 1) {
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
  # k
  if (!is.na(names(args["k"]))) {
    k <- eval(expression(k), parent.frame())
    if (!is.numeric(k) | length(k) != 1) {
      error_message(
        "k must be a number.",
        k,
        parent_function
      )
    }
    if (any(k < 1)) {
      error_message(
        "k must be positive",
        k,
        parent_function
      )
    }
    if (any(!is_integer_values(k))) {
      error_message(
        "k must be an integer",
        k,
        parent_function
      )
    }
  }
  # level
  if (!is.na(names(args["level"]))) {
    level <- eval(expression(level), parent.frame())
    if (!is.null(level)) {
      if (!is.numeric(level) | length(level) != 1) {
        error_message(
          "level must be a number.",
          level,
          parent_function
        )
      }
      if (any(level <= 0)) {
        error_message(
          "level must be positive.",
          level,
          parent_function
        )
      }
    }
  }
  # n
  if (!is.na(names(args["n"]))) {
    n <- eval(expression(n), parent.frame())
    if (!is.null(n)) {
      if (!is.numeric(n) | length(n) != 1) {
        error_message(
          "n must be a number.",
          n,
          parent_function
        )
      }
      if (any(n < 1)) {
        error_message(
          "n must be at least 1.",
          n,
          parent_function
        )
      }
    }
  }
  # n_simulations
  if (!is.na(names(args["n_simulations"]))) {
    n_simulations <- eval(expression(n_simulations), parent.frame())
    if (!is.numeric(n_simulations) | length(n_simulations) != 1) {
      error_message(
        "n_simulations must be a number.",
        n_simulations,
        parent_function
      )
    }
    if (any(n_simulations != 0 & n_simulations < 2)) {
      error_message(
        "n_simulations must be 0 or at least 2.",
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
  # orders
  if (!is.na(names(args["orders"]))) {
    orders <- eval(expression(orders), parent.frame())
    if (!is.null(orders)) {
      if (!is.numeric(orders) && !is.vector(orders)) {
        error_message(
          "orders must be a numeric vector",
          orders,
          parent_function
        )
      }
    }
  }
  # prob
  if (!is.na(names(args["prob"]))) {
    prob <- eval(expression(prob), parent.frame())
    if (!is.null(prob)) {
      if (!is.numeric(prob) && !is.vector(prob)) {
        error_message(
          "prob must be a numeric vector",
          prob,
          parent_function
        )
      }
      if (any(prob < 0)) {
        error_message(
          "probabilities must not be negative.",
          prob,
          parent_function
        )
      }
      if (abs(sum(prob) - 1) > length(prob) * .Machine$double.eps) {
        error_message(
          "probabilities must add up to 1.",
          prob,
          parent_function
        )
      }
    }
  }
  # probability_estimator is checked by match.arg()
  # prob_geom
  if (!is.na(names(args["prob_geom"]))) {
    prob_geom <- eval(expression(prob_geom), parent.frame())
    if (!is.null(prob_geom)) {
      if (!is.numeric(prob_geom) | length(prob_geom) != 1) {
        error_message(
          "prob_geom must be a number.",
          prob_geom,
          parent_function
        )
      }
      if (any(prob_geom < 0) | any(prob_geom > 1)) {
        error_message(
          "prob_geom must be between 0 and 1",
          prob_geom,
          parent_function
        )
      }
    }
  }
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
  # q_threshold
  if (!is.na(names(args["q_threshold"]))) {
    q_threshold <- eval(expression(q_threshold), parent.frame())
    if (!is.numeric(q_threshold) | length(q_threshold) != 1) {
      error_message(
        "q_threshold must be a number.",
        q_threshold,
        parent_function
      )
    }
    if (any(q_threshold < 1)) {
      error_message(
        "q_threshold must be positive",
        q_threshold,
        parent_function
      )
    }
    if (any(!is_integer_values(q_threshold))) {
      error_message(
        "q_threshold must be an integer",
        q_threshold,
        parent_function
      )
    }
  }
  # r
  if (!is.na(names(args["r"]))) {
    r <- eval(expression(r), parent.frame())
    if (!is.null(r)) {
      if (!is.numeric(r) && !is.vector(r)) {
        error_message(
          "r must be a numeric vector",
          r,
          parent_function
        )
      }
      if (length(r) < 2) {
        error_message(
          paste("r has length", length(r), "but it should be at least 2)"),
          r,
          parent_function
        )
      }
      if (r[1] != 0) {
        error_message(
          "The first r value must be 0",
          r,
          parent_function
        )
      }
      if (any(diff(r) <= 0)) {
        error_message(
          "Successive values of r must be increasing",
          r,
          parent_function
        )
      }
    }
  }
  # rate
  if (!is.na(names(args["rate"]))) {
    rate <- eval(expression(rate), parent.frame())
    if (!is.null(rate)) {
      if (!is.numeric(rate) | length(rate) != 1) {
        error_message(
          "rate must be a number.",
          rate,
          parent_function
        )
      }
      if (any(rate <= 0)) {
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
      if (!is.numeric(sample_coverage) | length(sample_coverage) != 1) {
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
  # sd_lnorm
  if (!is.na(names(args["sd_lnorm"]))) {
    sd_lnorm <- eval(expression(sd_lnorm), parent.frame())
    if (!is.null(sd_lnorm)) {
      if (!is.numeric(sd_lnorm) | length(sd_lnorm) != 1) {
        error_message(
          "sd_lnorm must be a number.",
          sd_lnorm,
          parent_function
        )
      }
      if (any(sd_lnorm <= 0)) {
        error_message(
          "sd_lnorm must be positive.",
          sd_lnorm,
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
  # size
  if (!is.na(names(args["size"]))) {
    size <- eval(expression(size), parent.frame())
    if (!is.null(size)) {
      if (!is.numeric(size) | length(size) != 1) {
        error_message(
          "size must be a number.",
          size,
          parent_function
        )
      }
      if (any(size <= 0)) {
        error_message(
          "size must be positive.",
          size,
          parent_function
        )
      }
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
  # species_number
  if (!is.na(names(args["species_number"]))) {
    species_number <- eval(expression(species_number), parent.frame())
    if (!is.null(species_number)) {
      if (!is.numeric(species_number) | length(species_number) != 1) {
        error_message(
          "species_number must be a number.",
          species_number,
          parent_function
        )
      }
      if (any(species_number < 1)) {
        error_message(
          "species_number must be at least 1.",
          species_number,
          parent_function
        )
      }
    }
  }
  # thomas_mu
  if (!is.na(names(args["thomas_mu"]))) {
    thomas_mu <- eval(expression(thomas_mu), parent.frame())
    if (!is.null(thomas_mu)) {
      if (!is.numeric(thomas_mu) | length(thomas_mu) != 1) {
        error_message(
          "thomas_mu must be a number.",
          thomas_mu,
          parent_function
        )
      }
      if (any(thomas_mu <= 0)) {
        error_message(
          "thomas_mu must be positive.",
          thomas_mu,
          parent_function
        )
      }
    }
  }
  # thomas_scale
  if (!is.na(names(args["thomas_scale"]))) {
    thomas_scale <- eval(expression(thomas_scale), parent.frame())
    if (!is.null(thomas_scale)) {
      if (!is.numeric(thomas_scale) | length(thomas_scale) != 1) {
        error_message(
          "thomas_scale must be a number.",
          thomas_scale,
          parent_function
        )
      }
      if (any(thomas_scale <= 0)) {
        error_message(
          "thomas_scale must be positive.",
          thomas_scale,
          parent_function
        )
      }
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
      if (!is.numeric(weights) && !is.vector(weights)) {
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
  # w_max
  if (!is.na(names(args["w_max"]))) {
    w_max <- eval(expression(w_max), parent.frame())
    if (!is.null(w_max)) {
      if (!is.numeric(w_max) | length(w_max) != 1) {
        error_message(
          "w_max must be a number.",
          w_max,
          parent_function
        )
      }
      if (any(w_max <= 0)) {
        error_message(
          "w_max must be positive.",
          w_max,
          parent_function
        )
      }
    }
  }
  # w_mean
  if (!is.na(names(args["w_mean"]))) {
    w_mean <- eval(expression(w_mean), parent.frame())
    if (!is.null(w_mean)) {
      if (!is.numeric(w_mean) | length(w_mean) != 1) {
        error_message(
          "w_mean must be a number.",
          w_mean,
          parent_function
        )
      }
      if (any(w_mean <= 0)) {
        error_message(
          "w_mean must be positive.",
          w_mean,
          parent_function
        )
      }
    }
  }
  # w_min
  if (!is.na(names(args["w_min"]))) {
    w_min <- eval(expression(w_min), parent.frame())
    if (!is.null(w_min)) {
      if (!is.numeric(w_min) | length(w_min) != 1) {
        error_message(
          "w_min must be a number.",
          w_min,
          parent_function
        )
      }
      if (any(w_min <= 0)) {
        error_message(
          "w_min must be positive.",
          w_min,
          parent_function
        )
      }
    }
  }
  # weibull_scale
  if (!is.na(names(args["weibull_scale"]))) {
    weibull_scale <- eval(expression(weibull_scale), parent.frame())
    if (!is.null(weibull_scale)) {
      if (!is.numeric(weibull_scale) | length(weibull_scale) != 1) {
        error_message(
          "weibull_scale must be a number.",
          weibull_scale,
          parent_function
        )
      }
      if (any(weibull_scale < 0)) {
        error_message(
          "weibull_scale must be positive.",
          weibull_scale,
          parent_function
        )
      }
    }
  }
  # weibull_shape
  if (!is.na(names(args["weibull_shape"]))) {
    weibull_shape <- eval(expression(weibull_shape), parent.frame())
    if (!is.null(weibull_shape)) {
      if (!is.numeric(weibull_shape) | length(weibull_shape) != 1) {
        error_message(
          "weibull_shape must be a number.",
          weibull_shape,
          parent_function
        )
      }
      if (any(weibull_shape <= 0)) {
        error_message(
          "weibull_shape must be strictly positive.",
          weibull_shape,
          parent_function
        )
      }
    }
  }
  # X
  if (!is.na(names(args["X"]))) {
    X <- eval(expression(X), parent.frame())
    if (!is.null(X)) {
      if (!inherits(X, "wmppp")) {
        error_message(
          "X must be an object of class 'wmppp'",
          X,
          parent_function
        )
      }
    }
  }
  # win
  if (!is.na(names(args["win"]))) {
    win <- eval(expression(win), parent.frame())
    if (!is.null(win)) {
      if (!inherits(X, "owin")) {
        error_message(
          "win must be an object of class 'owin'",
          win,
          parent_function
        )
      }
    }
  }

  # All tests passed.
  return(TRUE)
}


#' Chao's A
#'
#' Helper for Chao's estimators.
#'
#' A's formula \insertCite{@Chao2015@, eq. 6b}{divent}) depends on the presence
#' of singletons and doubletons.
#'
#' @param abd A vector of positive integers.
#'
#' @returns The value of A.
#' @references
#' \insertAllCited{}
#'
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
#' @returns A similarity matrix that corresponds to the species distribution.
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
#' @returns A tibble with the estimator used and the estimated entropy.
#' @noRd
#'
ent_gamma_tsallis <- function(
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
  if (q == 0 & estimator %in% c("jackknife", "iChao1", "Chao1", "rarefy", "naive")) {
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
  } else if (q == 1 & is.null(sample_coverage)) {
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
  } else if (q == 2 & is.null(sample_coverage)) {
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


#' Error message
#'
#' Utility for [check_divent_args]
#'
#' @param message The message to print.
#' @param argument The function argument that did not pass the tests.
#' @param parent_function The name of the function the argument was passed to.
#'
#' @returns Nothing. Used for side effects.
#' @noRd
#'
error_message <- function(message, argument, parent_function) {
  message(deparse(substitute(argument)), " cannot be:", appendLF = TRUE)
  message(utils::head(argument))
  message("In function ", parent_function, ": ", message)
  stop("Check the function arguments.", call. = FALSE)
}


#' Check Integers
#'
#' Check that the values of a vector are integer, whatever their type.
#'
#' @param x A numeric vector.
#'
#' @returns `TRUE` if values are integers.
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
#' @returns A list of matrices.
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
