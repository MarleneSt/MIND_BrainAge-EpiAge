###### install_packages.R: install the R packages used by the meta-analysis scripts ------
# Run once:   Rscript install_packages.R
# Versions used for the reported results are listed in README.md.

pkgs <- c("tidyverse", "metafor", "pheatmap", "Polychrome", "patchwork",
          "scales", "ggrepel", "corrplot", "cowplot")

to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install, repos = "https://cloud.r-project.org")

# check that all packages load
invisible(lapply(pkgs, library, character.only = TRUE))
cat("All packages installed and loaded.\n")
