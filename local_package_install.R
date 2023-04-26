## Installing R packages for the SCRGOT 2023 Coder Upgrade ##


# To install the R packages you need for the workshop, begin by opening RStudio
# and then copy the code below into the console.

# If you get a message about a package already being loaded, select "Yes" to
# restart R prior to installing the packages

packages_to_install_cran <-
    c("devtools",
      "markdown",
      "rmarkdown",
      "tidyverse",
      "Seurat",
      "qs",
      "patchwork",
      "msigdbr",
      "RColorBrewer",
      "ggrepel",
      "pheatmap",
      "knitr",
      "BiocManager",
      "SCpubr",
      "rslurm",
      "gridExtra")

install.packages(packages_to_install_cran, Ncpus = 4)

packages_to_install_bc <-
    c("fgsea",
      "harmony",
      "DESeq2",
      "SingleCellExperiment",
      "celldex",
      "SingleR",
      "limma",
      "ComplexHeatmap")

# You'll get a popup windows that says:
#   "Do you want to install from sources the packages which needs compilation?"
# Select "no"
# If it complains about out of date packages, select update none (type "n" and hit enter)
BiocManager::install(packages_to_install_bc, update = FALSE)


devtools::install_github("satijalab/seurat-data")
devtools::install_github("saeyslab/nichenetr")
devtools::install_github("rrrSingleCellUtils")
devtools::install_github("dynverse/dyno")


# If you get any messages about "had non-zero exit status", copy
# the error message and email it to me and I'll see if I can help.

