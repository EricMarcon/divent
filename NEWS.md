# divent 0.3-9

## Features

- Similarity and ordinariness of species.
- Similarity-based diversity: `ent_similarity()` and `div_similarity()`
- Rao's quadratic entropy: `ent_rao()`.

## TODO
- Check gamma = TRUE in div_ and ant_ functions, and in profile_ and accum_ functions.
- ent_phylo.species_distribution(gamma = TRUE) not implemented


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
