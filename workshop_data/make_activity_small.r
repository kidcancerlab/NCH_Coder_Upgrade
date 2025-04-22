library(Seurat)
library(DropletUtils) #enables us to write10xCounts
library(tidyverse)
library(Signac)
metadata <- read.table("Sessions/adv_thursday/activity_small_metadata.tsv")

activity_input_path <- "input/multiomics/downsampled_data/"

data_list <- qs::qread(paste0(activity_input_path, "gene_exp_peaks_list.qs"))

activity_ob <- CreateSeuratObject(counts = data_list[["Gene Expression"]],
                                  assay = "RNA",
                                  project = "10x_multiomics") %>%
    PercentageFeatureSet(pattern = "^MT",
                         col.name = "percent.mt_rna",
                         assay = "RNA")

activity_frag = paste0(activity_input_path, "new_fragments.tsv.gz")
activity_ob[["ATAC"]] <-
    CreateChromatinAssay(
        counts = data_list[["Peaks"]],
        sep = c(":", "-"),
        fragments = activity_frag,
        min.cells = 0
    )

activity_ob$bc <- colnames(activity_ob)

activity_small <- subset(activity_ob, bc %in% rownames(metadata))

activity_small@meta.data <- metadata
DefaultAssay(activity_small) <- "ATAC"
activity_small <- activity_small %>%
    RunTFIDF() %>%
    FindTopFeatures(min.cutoff = "q0") %>%
    RunSVD() %>%
    FindNeighbors(reduction = "lsi") %>%
    FindClusters(algorithm = 3) %>%
    RunUMAP(reduction = "lsi",
            dims = 2:30,
            reduction.name = "umap_atac")

DefaultAssay(activity_small) <- "RNA"
activity_small <- NormalizeData(activity_small) %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA() %>%
    RunUMAP(reduction.name = "umap_rna",
            dims = 1:10) %>%
    FindNeighbors(dims = 1:10) %>%
    FindClusters()

qs::qsave(activity_small, "input/multiomics/activity_data/activity_small.qs")