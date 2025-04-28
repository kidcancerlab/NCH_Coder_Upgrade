#!/usr/bin/env Rscript

# Notes on the RStudio server:
# Adam installed: quarto, libiconv, miniforge, python
# We needed to put LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib in environment
# We also increased /tmp to around 30G

# install with: dzdo Rscript pkg_installation.R

# Need to add /usr/local/lib to LD_LIBRARY_PATH so RCurl can find libiconv
# When running this script with dzdo the LD_LIBRARY_PATH isn't the same as my version
Sys.setenv(
    LD_LIBRARY_PATH = paste0(
        Sys.getenv("LD_LIBRARY_PATH"),
        ":/usr/local/lib"
    )
)

pkg_install_dir <- "/rstudio-workshop/apps/R/R-4.4.2_install/lib64/R/library/"

cran_dependencies <-
    c(
        "BiocManager",
        "pak",
        "fs"
    )

install.packages(
    cran_dependencies,
    Ncpus = 4,
    lib = pkg_install_dir,
    repos = "https://cran.rstudio.com"
)

bioc_dependencies <-
    c(
        "GenomeInfoDb",
        "GenomicRanges",
        "IRanges",
        "Rsamtools",
        "S4Vectors",
        "BiocGenerics",
        "limma",
        "ComplexHeatmap",
        "ensembldb",
        "SPOTlight",
        "biomaRt",
        "DropletUtils"
    )

needed_packages_cran <-
    c(
        "rstudioapi",
        "data.table",
        "DESeq2",
        "batchelor",
        "fgsea",
        "future",
        "ggrepel",
        "ggthemes",
        "harmony",
        "hdf5r",
        "knitr",
        "MAST",
        "msigdbr",
        "parallel",
        "patchwork",
        "pheatmap",
        "RColorBrewer",
        "Rcpp",
        "reticulate",
        "rmarkdown",
        "Seurat",
        "SeuratObject",
        "Signac",
        "testthat",
        "tidyverse",
        "terra",
        "anndata"
    )

pak::pkg_install(
    c(bioc_dependencies, needed_packages_cran),
    lib = pkg_install_dir
)

github_packages <-
    c(
        "saeyslab/nichenetr",
        "kidcancerlab/rrrSingleCellUtils",
        "NMikolajewicz/scMiko",
        "satijalab/seurat-data",
        "mojaveazure/seurat-disk",
        "satijalab/seurat-wrappers",
        "satijalab/azimuth",
        "drieslab/Giotto",
        "10xGenomics/loupeR",
        "jinworks/CellChat",
        "immunogenomics/crescendo",
        "dmcable/spacexr",
        "bnprks/BPCells/r"
    )

pak::pkg_install(
    github_packages,
    lib = pkg_install_dir
)

# Need to chmod the library folders afterwards as default seems to be 500 for pak
# ls -d /rstudio-workshop/apps/R/R-4.4.2_install/lib64/R/library/* | xargs dzdo chmod 555

warning(
    "Remember to run `ls -d /rstudio-workshop/apps/R/R-4.4.2_install/lib64/R/library/* | xargs dzdo chmod 555`"
)

# For rstudio server, modify the server's global settings
# [mvc002@rpl-rstudiows01 NCH_Coder_Upgrade]$ cat /etc/rstudio/rstudio-prefs.json
# {
#     "save_workspace": "never",
#     "load_workspace": false,
#     "knit_working_dir": "project",
#     "rmd_chunk_output_inline": false
# }

# Remember to chmod 644 on the file so user's instances can read it
