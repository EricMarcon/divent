#' Fit a distribution
#'
#' Fit a well-known distribution to a species distribution.
#'
#' [abundances] can be used to fit rank-abundance curves (RAC) of classical
#' distributions:
#'
#' - "lnorm" for log-normal \insertCite{Preston1948}{divent}.
#' - "lseries" for log-series \insertCite{Fisher1943}{divent}.
#' - "geom" for geometric \insertCite{Motomura1932}{divent}.
#' - "bstick" for broken stick \insertCite{MacArthur1957}{divent}.
#'    It has no parameter, so the maximum abundance is returned.
#'
#' @inheritParams check_divent_args
#' @param x An object
#' @param ... Unused.
#'
#' @returns A tibble with the sites and the estimated distribution parameters.
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' fit_rac(paracou_6_abd, distribution = "lnorm")
#'
#' @name fit_rac
NULL


#' @rdname fit_rac
#'
#' @export
fit_rac <- function(x, ...) {
  UseMethod("fit_rac")
}


#' @rdname fit_rac
#'
#'
#' @export
fit_rac.numeric <- function(
    x,
    distribution = c("lnorm", "lseries", "geom", "bstick"),
    ...,
    check_arguments = TRUE) {

  # Check arguments
  distribution <- match.arg(distribution)
  if (any(check_arguments)) {
    check_divent_args()
    if (any(x < 0)) stop("Species probabilities or abundances must be positive.")
  }

  # Eliminate zeros and sort.
  abd <- sort(x[x > 0], decreasing = TRUE)
  # Basic statistics
  s_obs <- length(abd)
  sample_size <- sum(abd)
  # Unique values
  nu <- unique(abd)

  if (distribution == "lnorm") {
    # Fit a lognormal distribution
    log_abd <- log(abd)
    mu <- mean(log_abd)
    sigma <- stats::sd(log_abd)
    rank <- s_obs * (1 - stats::pnorm(log(nu), mu, sigma))
    return(
      list(
        rac = tibble::tibble(
          rank = rank,
          abundance = nu
        ),
        parameters = tibble::tibble(
          mu = mu,
          sigma = sigma
        )
      )
    )
  }

  if (distribution == "lseries") {
    # Evaluate alpha
    alpha <- vegan::fisher.alpha(abd)
    # May (1975) Ecology and Evolution of Communities, Harvard University Press.
    sei <- function(t) exp(-t)/t
    rank <- vapply(
      nu,
      function(x) {
        n <- x * log(1 + alpha / sample_size)
        f <- stats::integrate(sei, n, Inf)
        fv <- f[["value"]]
        return(alpha * fv)
      },
      FUN.VALUE = 0
    )
    return(
      list(
        rac = tibble::tibble(
          rank = rank,
          abundance = nu
        ),
        parameters = tibble::tibble(
          alpha = alpha
        )
      )
    )
  }

  if (distribution == "geom") {
    # Fit a geometric distribution by a linear model
    log_abd <- log(abd)
    rank <- seq_len(s_obs)
    reg <- stats::lm(log_abd ~ rank)
    return(
      list(
        rac = tibble::tibble(
          rank = rank,
          abundance = exp(reg$coefficients[1] + reg$coefficients[2] * rank)
        ),
        parameters = tibble::tibble(
          prob = as.numeric(-reg$coefficients[2])
        )
      )
    )
  }

  if (distribution == "bstick") {
    # Fit a broken stick
    f1 <- sort(cumsum(1 / (s_obs:1)), decreasing = TRUE)
    nu <- sample_size * f1 / sum(f1)
    return(
      list(
        rac = tibble::tibble(
          rank = seq_len(s_obs),
          abundance = nu
        ),
        parameters = tibble::tibble(
          max = max(nu)
        )
      )
    )
  }
}


#' @rdname fit_rac
#'
#' @export
fit_rac.species_distribution <- function(
    x,
    distribution = c("lnorm", "lseries", "geom", "bstick"),
    ...,
    check_arguments = TRUE) {

  # Check arguments
  distribution <- match.arg(distribution)
  if (any(check_arguments)) {
    check_divent_args()
    if (any(x < 0)) stop("Species probabilities or abundances must be positive.")
  }

  # Apply probabilities.numeric() to each site
  rac_list <- apply(
    # Eliminate site and weight columns
    x[, !colnames(x) %in% non_species_columns],
    # Apply to each row
    MARGIN = 1,
    FUN = fit_rac.numeric,
    # Arguments
    distribution = distribution,
    check_arguments = FALSE
  )

  # Bind the rows of the parameters tibble
  racs <- dplyr::bind_rows(lapply(rac_list, function(x) x$parameters))
  # Restore the site names
  if ("site" %in% colnames(x)) {
    racs <- dplyr::bind_cols(site = x$site, racs)
  }
  return(racs)
}
