cran_packages <- c(
  "dplyr", "tidyr", "data.table", "cowplot", "igraph", "purrr", "BiocManager"
)

bioc_packages <- c(
  "GenomicRanges", "IRanges"
)

install_if_missing <- function(pkgs, installer) {
  missing <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(missing) > 0) {
    message("Installing missing packages: ", paste(missing, collapse = ", "))
    installer(missing)
  }
}

install_if_missing(cran_packages, install.packages)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install_if_missing(bioc_packages, BiocManager::install)

message("Environment setup complete.")

cran_packages <- c(
  "dplyr", "tidyr", "data.table", "cowplot", "igraph", "purrr"
)

bioc_packages <- c(
  "GenomicRanges", "IRanges"
)

all_packages <- c(cran_packages, bioc_packages)

invisible(lapply(all_packages, function(pkg) {
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE)
  )
}))
