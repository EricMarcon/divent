# divent 0.5-1

## Features

- Generalized Simpson's entropy: `ent_gen_simpson()` and `div_gen_simpson()`
- Faith's PD: `div_pd()`
- Allen at al's phylogenetic entropy: `ent_allen()`
- Spatially explicit Simpson's entropy: `ent_sp_simpson()` and `ent_sp_simpsonEnvelope()`
- Spatially explicit random communities: `rspcommunity()`
- `species_distribution` methods for `wmppp` and `character` objects.
- values of abundances are no longer limited to `.Machine$integer.max`.

## External changes

- Replaced `Geom*$default_aes` by their values for compatibility with ggplot2 3.6.0 (PR #2 by @teunbrand)


# divent 0.4-4

## Features

- Hurlbert's diversity: `ent_hurlbert()` and `div_hurlbert()`.


# divent 0.3-16

## Features

- Similarity and ordinariness of species.
- Similarity-based diversity: `ent_similarity()` and `div_similarity()`
- Rao's quadratic entropy: `ent_rao()`.
- paracou_6_wmppp dataset.


# divent 0.2-5

## Features

- Data formats: 
    - phylo_divent (phylogenetic tree): `as_phylo_divent()`
- Phylogenetic diversity: `ent_phylo()` and `div_phylo()`
- species names must be valid names


# divent 0.1-24

- First Version.

## Features

- Data formats: 
    - species distribution: `species_distribution()`
    - metacommunity: `metacommunity()`
- Sample coverage: `coverage()`
- Richness: `div_richness()`
- Shannon's, Simpson's and Tsallis's entropies: `ent_shannon()`, `ent_simpson()`, `ent_tsallis()`
- Hill numbers: `div_hill()`
- Diversity accumulation: `accum_hill()`
- Diversity profiles: `profile_hill()`
- Diversity partitioning: `profile_hill()`
