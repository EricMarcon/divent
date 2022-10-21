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
#' @references
#' \insertAllCited{}
NULL

# Names of variables:
# abd: a numeric vector of abundances
# prob: a numeric vector of probabilities
# abundances / probabilities: an object of class abundances / probabilities
# s_0, s_1,  ...: species observed 0, 1, ... times
# sample_size: number of observed individuals
# species_number: number of observed species
# coverage_deficit_2: coverage deficit of order 2 of the sample
# sample_coverage: coverage of order 1


#  Utilities ----
#' check_divent_args
#'
#' Checks the arguments of a function of the package divent
#'
#' The function compares the arguments passed to its parent function to the type 
#' they should be and performs some extra tests, *e.g.* probabilities must be positive and sum to 1. 
#' It stops if an argument is not correct.
#'
#' @return Returns `TRUE` or stops if a problem is detected.
#' 
#' @export
#'
#' @keywords internal
#' 
check_divent_args <- function() {

  # Verify that the package is attached
  if (! "divent" %in% .packages()) {
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
  
  # q
  if (!is.na(names(args["q"]))) {
    q <- eval(expression(q), parent.frame())
    if (!is.numeric(q) | length(q) != 1)
      error_message("q must be a number.", q, parent_function)
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
