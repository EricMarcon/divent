# Changelog

## divent 0.5-3

CRAN release: 2025-08-27

### Bug correction

- [`species_distribution()`](https://ericmarcon.github.io/divent/reference/species_distribution.md)
  raised an error when input was a vector of factors.

## divent 0.5-2

CRAN release: 2025-02-09

### Features

- Generalized Simpson’s entropy:
  [`ent_gen_simpson()`](https://ericmarcon.github.io/divent/reference/ent_gen_simpson.md)
  and
  [`div_gen_simpson()`](https://ericmarcon.github.io/divent/reference/div_gen_simpson.md)
- Faith’s PD:
  [`div_pd()`](https://ericmarcon.github.io/divent/reference/div_pd.md)
- Allen at al’s phylogenetic entropy:
  [`ent_allen()`](https://ericmarcon.github.io/divent/reference/ent_allen.md)
- Spatially explicit Simpson’s entropy:
  [`ent_sp_simpson()`](https://ericmarcon.github.io/divent/reference/ent_sp_simpson.md)
  and
  [`ent_sp_simpsonEnvelope()`](https://ericmarcon.github.io/divent/reference/ent_sp_simpson.md)
- Spatially explicit random communities:
  [`rspcommunity()`](https://ericmarcon.github.io/divent/reference/rcommunity.md)
- `species_distribution` methods for `wmppp` and `character` objects.
- values of abundances are no longer limited to `.Machine$integer.max`.

### External changes

- Replaced `Geom*$default_aes` by their values for compatibility with
  ggplot2 3.6.0 (PR [\#2](https://github.com/EricMarcon/divent/issues/2)
  by [@teunbrand](https://github.com/teunbrand))

## divent 0.4-4

CRAN release: 2024-11-06

### Features

- Hurlbert’s diversity:
  [`ent_hurlbert()`](https://ericmarcon.github.io/divent/reference/ent_hurlbert.md)
  and
  [`div_hurlbert()`](https://ericmarcon.github.io/divent/reference/div_hurlbert.md).

## divent 0.3-16

### Features

- Similarity and ordinariness of species.
- Similarity-based diversity:
  [`ent_similarity()`](https://ericmarcon.github.io/divent/reference/ent_similarity.md)
  and
  [`div_similarity()`](https://ericmarcon.github.io/divent/reference/div_similarity.md)
- Rao’s quadratic entropy:
  [`ent_rao()`](https://ericmarcon.github.io/divent/reference/ent_rao.md).
- paracou_6_wmppp dataset.

## divent 0.2-5

### Features

- Data formats:
  - phylo_divent (phylogenetic tree):
    [`as_phylo_divent()`](https://ericmarcon.github.io/divent/reference/phylo_divent.md)
- Phylogenetic diversity:
  [`ent_phylo()`](https://ericmarcon.github.io/divent/reference/ent_phylo.md)
  and
  [`div_phylo()`](https://ericmarcon.github.io/divent/reference/div_phylo.md)
- species names must be valid names

## divent 0.1-24

- First Version.

### Features

- Data formats:
  - species distribution:
    [`species_distribution()`](https://ericmarcon.github.io/divent/reference/species_distribution.md)
  - metacommunity:
    [`metacommunity()`](https://ericmarcon.github.io/divent/reference/metacommunity.md)
- Sample coverage:
  [`coverage()`](https://ericmarcon.github.io/divent/reference/coverage.md)
- Richness:
  [`div_richness()`](https://ericmarcon.github.io/divent/reference/div_richness.md)
- Shannon’s, Simpson’s and Tsallis’s entropies:
  [`ent_shannon()`](https://ericmarcon.github.io/divent/reference/ent_shannon.md),
  [`ent_simpson()`](https://ericmarcon.github.io/divent/reference/ent_simpson.md),
  [`ent_tsallis()`](https://ericmarcon.github.io/divent/reference/ent_tsallis.md)
- Hill numbers:
  [`div_hill()`](https://ericmarcon.github.io/divent/reference/div_hill.md)
- Diversity accumulation:
  [`accum_hill()`](https://ericmarcon.github.io/divent/reference/accum_hill.md)
- Diversity profiles:
  [`profile_hill()`](https://ericmarcon.github.io/divent/reference/profile_hill.md)
- Diversity partitioning:
  [`profile_hill()`](https://ericmarcon.github.io/divent/reference/profile_hill.md)
