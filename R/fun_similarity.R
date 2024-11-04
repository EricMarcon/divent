#' Functional similarity
#' 
#' Transform a distance matrix into a similarity matrix 
#' \insertCite{Leinster2012}{divent}.
#' Similarity between two species is defined either by a negative exponential 
#' function of their distance or by the complement to 1 of their normalized 
#' distance (such that the most distant species are 1 apart).
#'
#' @inheritParams check_divent_args
#' @param exponential If `TRUE`, similarity is \eqn{e^{-r \delta}}, 
#' where \eqn{r} is argument `rate`.
#' If `FALSE`, it is \eqn{1 - \delta / \max(\delta)}.
#'
#' @returns A similarity matrix.
#' @export
#'
#' @references
#' \insertAllCited{}
#' 
#' @examples
#' # Similarity between Paracou 6 species
#' hist(fun_similarity(paracou_6_fundist))
#' 
fun_similarity <- function (
    distances,
    exponential = TRUE,
    rate = 1,
    check_arguments = TRUE) {
  
  if (any(check_arguments)) {
    check_divent_args()
    # Names
    if (xor(is.null(colnames(distances)), is.null(rownames(distances)))) {
      stop("Row and column names must be identical in 'distances'.")
    }
  }
  
  if (inherits(distances, "dist")) {
    # dist objects are supported but the remainder assumes a matrix
    distances <- as.matrix(distances)
  }

  if (is.null(colnames(distances))) {
    ### Add default species names such as sp_1
    colnames(x) <- paste(
      "sp", 
      formatC(seq_len(ncol(x)), width = ceiling(log10(ncol(x))), flag = "0"),
      sep = "_"
    )
    rownames(x) <- colnames(x)
  } else if (!identical(colnames(distances), rownames(distances))) {
    stop("Row and column names must be identical in 'distances'.")
  }
  
  if (exponential) {
    the_similarity <- exp(-rate * distances)
  } else {
    the_similarity <- 1 - distances / max(distances)
  }
  
  class(the_similarity) <- c("similarity", class(the_similarity))
  return(the_similarity)
}
