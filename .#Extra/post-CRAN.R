# To be run after each CRAN update

# Install packages ----
packages_required <- c("codemetar", "cffr")
# Already installed?
packages_installed <- vapply(
  packages_required,
  FUN = requireNamespace,
  FUN.VALUE = TRUE,
  quietly = TRUE
)
# Install missing packages
install.packages(
  names(packages_installed)[!packages_installed],
  repos = "https://cran.rstudio.com/",
  quiet = TRUE
)

# Update codemeta.json ----
codemetar::write_codemeta(".")

# Update CITATION.cff ----
cffr::cff_write()
