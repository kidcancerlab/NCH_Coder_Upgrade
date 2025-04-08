#exploring activity dataset
# input/multiomics/Duan_029_20088-1_input_lib_hg38_only
library(tidyverse)
library(Seurat)
library(rrrSingleCellUtils)
library(Signac)

activity_input_path <-
    "input/multiomics/Duan_029_20088-1_input_lib_hg38only/outs/"

h5_list <- Read10X_h5(paste0(activity_input_path,
                             "filtered_feature_bc_matrix.h5"))
set.seed(17)
down_index <- sample(1:ncol(h5_list[[1]]), 2000) #noling
names(h5_list)
activity_ob <-
    CreateSeuratObject(counts = h5_list[["Gene Expression"]][ , down_index],
                       assay = "RNA",
                       project = "10x_multiomics") %>%
    PercentageFeatureSet(pattern = "^MT",
                         col.name = "percent.mt_rna",
                         assay = "RNA")

#Need to downsample this
#think I'll downsample to ~2k, so I still have enough for QC
activity_ob <- process_seurat(activity_ob)
DimPlot(activity_ob)

#going to add atac assay
head(rownames(h5_list[["Peaks"]]))
activity_frag <- paste0(activity_input_path, "atac_fragments.tsv.gz")
activity_ob[["ATAC"]] <-
    CreateChromatinAssay(counts = h5_list[["Peaks"]][ , down_index],
                         sep = c(":", "-"),
                         fragments = activity_frag,
                         min.cells = 0)

#Add gene annotation information
DefaultAssay(activity_ob) <- "ATAC"
annotations <-
    GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
seqlevels(annotations) <- paste0("chr", seqlevels(annotations))
genome(annotations) <- "hg38"

Annotation(activity_ob) <- annotations

#maybe scrap activity 1 and just provide object up till this point? and maybe add in the FRiP since we don't have summary info, or show other way to calculate frip
#Add nucleosome signal ~3 seconds
start <- Sys.time()
activity_ob <- NucleosomeSignal(activity_ob, assay = "ATAC")
end <- Sys.time()
end - start

#TSS enrichment ~3.5 minutes
start <- Sys.time()
activity_ob <- TSSEnrichment(activity_ob, assay = "ATAC", fast = FALSE)
Sys.time() - start

#FRiP
## counting fragments 1.5 minutes
start <- Sys.time()
total_frag_df <- CountFragments(activity_frag) %>%
    dplyr::filter(CB %in% colnames(activity_ob)) %>%
    arrange(match(CB, colnames(activity_ob)))
Sys.time() - start

activity_ob@meta.data <- activity_ob@meta.data %>%
    dplyr::mutate(total_frag = total_frag_df$frequency_count*2,
                  mononucleosomal = total_frag_df$mononucleosomal,
                  nucleosome_free = total_frag_df$nucleosome_free)
activity_ob <- FRiP(activity_ob,
                    assay = "ATAC",
                    total.fragments = "total_frag",
                    col.name = "FRiP")

#plot qc metrics
Idents(activity_ob) <- NULL
VlnPlot(activity_ob,
        features = c("nCount_RNA",
                     "percent.mt_rna",
                     "nucleosome_signal",
                     "tss.enrichment",
                     "FRiP"),
        group.by = "orig.ident",
        ncol = 4)
