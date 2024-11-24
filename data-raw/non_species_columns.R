## code to prepare `non_species_columns` variable

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


# Save
usethis::use_data(non_species_columns, overwrite = TRUE)
