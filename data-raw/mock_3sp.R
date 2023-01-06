## code to prepare `mock_3sp` dataset

# Distance matrix
mock_3sp <- matrix(c(0, 1, 2, 1, 0, 2, 2, 2, 0), nrow = 3, byrow = TRUE)
row.names(mock_3sp) <- colnames(mock_3sp) <- c("A", "B", "C")
mock_3sp

# dist object
mock_3sp_dist <- as.dist(mock_3sp)

# phylo tree
mock_3sp_tree <- ape::as.phylo(hclust(mock_3sp_dist, method = "average"))
plot(mock_3sp_tree); axis(1)

# Abundance vector
mock_3sp_abd <- 3:1
names(mock_3sp_abd) <- labels(mock_3sp_dist)

# Save
usethis::use_data(mock_3sp_abd, overwrite = TRUE)
usethis::use_data(mock_3sp_dist, overwrite = TRUE)
usethis::use_data(mock_3sp_tree, overwrite = TRUE)
