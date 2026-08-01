# Todo list for divent as of version 0.5.4

## Documentation

- Detail the estimators in `div_similarity()` and `ent_similarity()`.

## Code

### `div_part()`

1. Support a vector of orders.

2. Support phylo and similarity diversity: 
  - add arguments tree and Z (entropart style)
  - or add function `div_part_phylo()` and `div_part_similarity()` (better, but div_part should be renamed `div_part_hill()`).
  - or combine both: `div_part()` as a wrapper for the three functions.

### Plot diversity partitioning profiles

- Add an `autoplot()` able similar to that of entropart.

- Add `profile_beta_*()` functions. Accept an object created by `div_part()` as input.

### Spatial Accumulation of Diversity

- Rename that Neighborhood diversity, cite Rao's tribute book chapter.
Rename `accum_sp_*()` functions `nbd_*()` for clarity (avoid confusion with DACs `accum_*()`.

- add `nbd_phylo()` and `nbd_similarity()` functions.

- parallelize the computing of the entropy of each neigborhood.

### ISARs

- add `isar()` function accepting individual neighborhood diversity as input.

### Neighborhood entropy

Complete `ent_sp_simpson()` with `ent_sp_rao()` and link to package ads for the standardized versions.

