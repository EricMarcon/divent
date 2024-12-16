#' Spatially Explicit Simpson's Entropy
#'
#' Simpson's entropy of the neighborhood of individuals,
#' up to a distance \insertCite{Shimatani2001}{divent}.
#'
#' @inheritParams check_divent_args
#' @param correction the edge-effect correction to apply when estimating
#' the number of neighbors or the *K* function with [spatstat.explore::Kest].
#' Default is "isotropic".
#'
#' @name ent_sp_simpson
#' @references
#' \insertAllCited{}
NULL


#' @rdname ent_sp_simpson
#' @returns `ent_sp_simpson` returns an object of class `fv`,
#' see [spatstat.explore::fv.object].
#' There are methods to print and plot this class.
#' It contains the value of the spatially explicit Simpson's entropy
#' for each distance in `r`.
#' @export
#'
#' @examples
#' # Generate a random community
#' X <- rspcommunity(1, size = 1000, species_number = 3)
#' # Calculate the entropy and plot it
#' autoplot(ent_sp_simpson(X))
#'
ent_sp_simpson <- function(
    X,
    r = NULL,
    correction = c("isotropic", "translate", "none"),
    check_arguments = TRUE) {

  # Check arguments
  correction <- match.arg(correction)
  if (any(check_arguments)) {
    check_divent_args()
  }

  # Summary
  abd <- tapply(
    spatstat.geom::marks(X)$PointType,
    spatstat.geom::marks(X)$PointType,
    length)
  abd <- abd[!is.na(abd)]
  sample_size <- sum(abd)
  prob <- abd / sample_size

  # r
  if (is.null(r)) {
    r_max <- spatstat.geom::diameter(X$window)
    # Default values, as in dbmss
    r <- r_max *
      c(
        0,
        1:20,
        seq(22, 40, 2),
        seq(45, 100, 5),
        seq(110, 200, 10),
        seq(220, 400, 20)
      ) / 800
  }
  # K all points
  K_all <- correction_fv(
    spatstat.explore::Kest(
      X,
      r = r,
      correction = correction
    ),
    correction
  )

  # The point pattern is split into a list of ppp for each mark
  ppp.list <- split(X, as.factor(spatstat.geom::marks(X)$PointType))
  # K for each ppp
  K.list <- lapply(
    ppp.list,
    FUN = spatstat.explore::Kest,
    r = r,
    correction = correction
  )
  # Extract the values into a dataframe, each column is a PointType
  K_PointType <- as.data.frame(
    lapply(
      K.list,
      FUN = correction_fv,
      correction = correction
    )
  )
  # K_PointType is NA for species with a a single point. Should be 0
  K_PointType[is.na(K_PointType)] <- 0

  # Result: Shilmatani's function of r
  Shi_r <- (1 - rowSums(
    (K_PointType * rep(abd * (abd - 1), each = dim(K_PointType)[1])) /
      (K_all * sample_size * (sample_size - 1)))
  ) * (sample_size - 1) / sample_size

  # Build a dataframe with r, theoretical value and S(r)
  Shi.df <- data.frame(
    r = r,
    Simpson = ent_simpson.numeric(abd, as_numeric = TRUE),
    S_r = Shi_r
  )

  # Return the values of Shimatani(r), aka Simpson(r)
  labl <- c("r", "hat(%s)", "hat(%s)(r)")
  desc <- c("Distance argument r", "Asymptotic %s", "Estimated %s")
  the_entropy <- spatstat.explore::fv(
    Shi.df,
    argu = "r",
    ylab = quote(Shimatani(r)),
    valu = "S_r",
    fmla = ". ~ r",
    alim = c(0, max(r)),
    labl = labl,
    desc = desc,
    unitname = X$window$unit,
    fname = "Simpson's Entropy")
  spatstat.explore::fvnames(the_entropy, ".") <- colnames(Shi.df)[-1]
  return(the_entropy)
}



#' @rdname ent_sp_simpson
#' @param h0 A string describing the null hypothesis to simulate.
#' The null hypothesis may be "RandomPosition": points are drawn in a Poisson process (default)
#' or "RandomLabeling": randomizes point types, keeping locations unchanged.
#'
#' @return `ent_sp_simpsonEnvelope` returns an envelope object [spatstat.explore::envelope].
#' There are methods to print and plot this class.
#' It contains the observed value of the function,
#' its average simulated value and the confidence envelope.
#' @export
#'
#' @examples
#' # Generate a random community
#' X <- rspcommunity(1, size = 1000, species_number = 3)
#' # Calculate the entropy and plot it
#' autoplot(ent_sp_simpsonEnvelope(X, n_simulations = 10))
#'
ent_sp_simpsonEnvelope <- function(
    X,
    r = NULL,
    n_simulations = 100,
    alpha = 0.05,
    correction = c("isotropic", "translate", "none"),
    h0 = c("RandomPosition", "RandomLabeling"),
    global = FALSE,
    check_arguments = TRUE) {

  # Check arguments
  correction <- match.arg(correction)
  h0 <- match.arg(h0)
  if (any(check_arguments)) {
    check_divent_args()
  }

  # Choose the null hypothesis
  X_sim <- switch(
    h0,
    RandomPosition = expression(dbmss::rRandomPositionK(X, CheckArguments = FALSE)),
    RandomLabeling = expression(dbmss::rRandomLabeling(X, CheckArguments = FALSE))
  )
  if (is.null(X_sim)) {
    stop(paste("The null hypothesis", sQuote(h0), "has not been recognized."))
  }
  # local envelope, keep extreme values for lo and hi (nrank=1)
  the_envelope <- spatstat.explore::envelope(
    X,
    fun = ent_sp_simpson,
    nsim = n_simulations,
    nrank = 1,
    r = r,
    correction = correction,
    check_arguments = FALSE,
    simulate = X_sim,
    savefuns = TRUE
  )
  attr(the_envelope, "einfo")$H0 <- switch(
    h0,
    RandomPosition = "Random Position",
    RandomLabeling = "Random Labeling"
  )
  # Calculate confidence intervals
  the_envelope <- dbmss::FillEnvelope(the_envelope, Alpha = alpha, Global = global)
  # Return the envelope
  return(the_envelope)
}


#' Extract a column from an fv object
#' according to an edge-effect correction
#'
#' @param fv the function value object, see [spatstat.explore::fv.object].
#' @param correction the edge-effect correction:
#' "isotropic", "translate" or "none"
#'
#' @returns a vector with the function values
#' @noRd
#'
correction_fv <- function(fv, correction) {
  switch(
    correction,
    "isotropic" = fv$iso,
    "translate" = fv$trans,
    "none" = fv$un
  )
}
