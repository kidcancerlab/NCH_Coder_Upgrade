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

## Libraries needed to be installed:
# libuv
# xml2
# libgsl
# libgdal
# libgeos
# libproj
# libtbb
# libnetcdf
# hdf5
# fftw3
# freetype2
# fontconfig
# CMake
# NLopt
# fribidi

### libjpg
### libwebp  dzdo yum install libjpeg-devel libwebp-devel
### cairo    dzdo yum install cairo-devel

# Make sure this matches .libPaths()
# Also double check the chmod at the end
pkg_install_dir <- "/rstudio-workshop/apps/R/R-4.5.3/lib64/R/library"

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
    "DropletUtils",
    "multtest",
    "metap",
    "biovizBase",
    "multtest",
    "BSgenome.Mmusculus.UCSC.mm10",
    "BSgenome.Mmusculus.UCSC.mm39"
  )

needed_packages_cran <-
  c(
    "devtools",
    "rstudioapi",
    "data.table",
    "DESeq2",
    "batchelor",
    "future",
    "ggrepel",
    "ggthemes",
    "harmony",
    "hdf5r",
    "knitr",
    "MAST",
    "patchwork",
    "pheatmap",
    "RColorBrewer",
    "Rcpp",
    "reticulate",
    "rmarkdown",
    "Seurat",
    "SeuratObject",
    "testthat",
    "tidyverse",
    "terra",
    "anndata",
    "scatterpie",
    "qs2",
    "scrapper"
  )

#Doing this one at a time so I can see which one fails
for (this_pkg in c(bioc_dependencies, needed_packages_cran)) {
  message("Now installing: ", this_pkg)

  pak::pkg_install(
    this_pkg,
    lib = pkg_install_dir
  )
}

# I need to install a specific version of signac to get azimuth to install right
remotes::install_github("stuart-lab/signac", ref = "1.16.0")

github_packages <-
  c(
    "igordot/msigdbr",
    "alserglab/fgsea",
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
    "JEFworks-Lab/MERINGUE",
    "KlugerLab/ALRA",
    "jinworks/SpatialCellChat",
    "immunogenomics/crescendo",
    "dmcable/spacexr",
    "MVesuviusC/ramTrackR"
  )

for (this_pkg in github_packages) {
  message("Now installing: ", this_pkg)

  pak::pkg_install(
    this_pkg,
    lib = pkg_install_dir
  )
}

# Need to get a specific version of harmony to get around an issue (https://github.com/immunogenomics/harmony/issues/284)
devtools::install_github("pati-ni/harmony", ref="lapack-linux-fix", force=TRUE)

# This failed when installed with pak
# Error was
# ! error in pak subprocess
# Caused by error in `stop_task_package_uncompress(state, worker)`:
#   ! Failed to uncompress BPCells from
# /tmp/RtmpMXar8O/file659fc295acc10/src/contrib/BPCells_0.3.1_bde5737.tar.gz-t.
# Similar issue: https://github.com/bnprks/BPCells/issues/326
devtools::install_github("bnprks/BPCells/r", lib = pkg_install_dir)


# motifmatcher was complaining about the default C++ compiler being C++11 instead
# of C++14. I had to force R to use 14
# I did this at the end so it doesn't mess up anything else
tmp <- tempfile("Makevars")
writeLines(
  'CXX11 = g++
CXX11STD = -std=gnu++14
CXX14 = g++
CXX14STD = -std=gnu++14',
  con = tmp
)
Sys.setenv(R_MAKEVARS_USER = tmp)

BiocManager::install("motifmatchr")


# Need to chmod the library folders afterwards as default seems to be 500 for pak
# ls -d /rstudio-workshop/apps/R/R-4.5.3/lib64/R/library/* | xargs dzdo chmod 555

warning(
  "Remember to run `ls -d /rstudio-workshop/apps/R/R-4.5.3/lib64/R/library/* | xargs dzdo chmod 555`"
)

# For rstudio server, modify the server's global settings
# [mvc002@rpl-rstudiows01 NCH_Coder_Upgrade]$ cat /etc/rstudio/rstudio-prefs.json
# {
#     "save_workspace": "never",
#     "load_workspace": false,
#     "knit_working_dir": "project",
#     "rmd_chunk_output_inline": false,
#     "graphics_backend": "ragg"
# }

# Remember to chmod 644 on the file so user's instances can read it
