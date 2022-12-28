# divent 0.2-3

## Features

- Data formats: 
    - phylo_divent (phylogenetic tree): `as_phylo_divent()`
- Phylogenetic entropy: `ent_phylo()`


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
- Diversity accumulation: `div_accum()`
- Diversity profiles: `div_profile()`
- Diversity partitioning: `div_profile()`

## TODO
- `div_profile` requires integer abundances for its confidence envelope but `gamma` argument is ignored.
