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

activity_ob <-
    CreateSeuratObject(counts = h5_list[["Gene Expression"]],
                       assay = "RNA",
                       project = "10x_multiomics") %>%
    PercentageFeatureSet(pattern = "^MT",
                         col.name = "percent.mt_rna",
                         assay = "RNA")

#going to add atac assay
head(rownames(h5_list[["Peaks"]]))
activity_frag <- paste0(activity_input_path, "atac_fragments.tsv.gz")
activity_ob[["ATAC"]] <-
    CreateChromatinAssay(counts = h5_list[["Peaks"]],
                         sep = c(":", "-"),
                         fragments = activity_frag,
                         min.cells = 0)

# Need to downsample this
set.seed(17)
bcs_keep <- sample(colnames(activity_ob), size = 2000)
activity_ob$keep_cell <- colnames(activity_ob) %in% bcs_keep
activity_ob <- subset(activity_ob, keep_cell)

qs::qsave(activity_ob, "testing_folder/downsampled_ob.qs")

#think I'll downsample to ~2k, so I still have enough for QC
DimPlot(activity_ob)

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

qs::qsave(activity_ob, "testing_folder/ob_w_qc_metrics.qs")

#plot qc metrics
VlnPlot(activity_ob,
        features = c("nCount_RNA",
                     "percent.mt_rna",
                     "nucleosome_signal",
                     "TSS.enrichment",
                     "FRiP"),
        group.by = "orig.ident",
        ncol = 5)

cutoffs <- tribble(~feature, ~min_val, ~max_val,
                   "nCount_RNA", 1000, 20000,
                   "percent.mt_rna", 0, 20,
                   "nucleosome_signal", 0, 2,
                   "TSS.enrichment", 1, 50,
                   "FRiP", 0.25, 1)

feature_hist(activity_ob, features = cutoffs$feature, cutoff_table = cutoffs)

activity_ob <- subset(activity_ob,
                      nCount_RNA %in% c(1000:20000) &
                      percent.mt_rna < 20 &
                      nucleosome_signal < 2 &
                      TSS.enrichment > 1 &
                      FRiP > 0.25)

DefaultAssay(activity_ob) <- "RNA"
activity_ob <- NormalizeData(activity_ob) %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA() %>%
    RunUMAP(reduction.name = "umap_rna",
            dims = 1:10) %>%
    FindNeighbors(dims = 1:10) %>%
    FindClusters()

activity_ob$RNA_cluster <- Idents(activity_ob)

DefaultAssay(activity_ob) <- "ATAC"
activity_ob <- RunTFIDF(activity_ob) %>%
    FindTopFeatures(min.cutoff = "q0") %>%
    RunSVD() %>%
    FindNeighbors(reduction = "lsi") %>%
    FindClusters(algorithm = 3) %>%
    RunUMAP(reduction = "lsi",
            dims = 2:30,
            reduction.name = "umap_atac")

activity_ob$ATAC_cluster <- Idents(activity_ob)

qs::qsave(activity_ob, "testing_folder/qced_activity_ob.qs")

r_dim_plot(activity_ob, reduction = "umap_rna", group.by = "RNA_cluster") |
r_dim_plot(activity_ob, reduction = "umap_atac", group.by = "ATAC_cluster")


#Last activity
#make small object with just clusters 0 and 1
activity_small <- subset(activity_ob, ATAC_cluster %in% c(0, 1)) %>%
    RunTFIDF() %>%
    FindTopFeatures(min.cutoff = "q0") %>%
    RunSVD() %>%
    FindNeighbors(reduction = "lsi") %>%
    FindClusters(algorithm = 3) %>%
    RunUMAP(reduction = "lsi",
            dims = 2:30,
            reduction.name = "umap_atac")

qs::qsave(activity_small, "testing_folder/activity_small.qs")

#find differentially accessible peaks
Idents(activity_small) <- activity_small$ATAC_cluster

r_dim_plot(activity_small, reduction = "umap_rna", group.by = "RNA_cluster") |
r_dim_plot(activity_small, reduction = "umap_atac", group.by = "ATAC_cluster")

diff_peaks <- FindAllMarkers(activity_small,
                          assay = "ATAC",
                          min.pct = .1) %>%
    subset(p_val_adj < 0.05)

top_5_each <- group_by(diff_peaks, cluster) %>%
    arrange(desc(avg_log2FC)) %>%
    slice_head(n = 5)

#plot peaks as atac feature plot
r_feature_plot(activity_ob, features = top_5_each$gene)

#Link peaks to genes
#First I'm going to find differentially expressed genes between ATAC clusters, bc we're only interested in genes whose expression is changing between states
DefaultAssay(activity_small) <- "RNA"
activity_small <- process_seurat(activity_small)
activity_small$RNA_cluster <- activity_small$seurat_clusters
Idents(activity_small) <- activity_small$RNA_cluster
degs <- FindAllMarkers(activity_small,
                       assay = "RNA",
                       only.pos = TRUE,
                       min.pct = 0.1) %>%
    subset(p_val_adj < 0.05 & avg_log2FC > 1)
dim(degs)

#link peaks to genes; this returns a seurat object
#have to run region stats first
DefaultAssay(activity_small) <- "ATAC"
start <- Sys.time()
activity_small <-
    RegionStats(activity_small,
                genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38) #nolint
Sys.time() - start

#Link peaks takes about 7 minutes
start <- Sys.time()
activity_small <- LinkPeaks(object = activity_small,
                            peak.assay = "ATAC",
                            expression.assay = "RNA",
                            genes.use = rownames(degs)[1:200])
Sys.time() - start
qs::qsave(activity_small, "testing_folder/activity_small.qs")

#look at peaks
links <- Links(activity_small) %>%
    as.data.frame() %>%
    as_tibble() %>%
    arrange(pvalue)
head(links)

#lets see if any of our top 5 peaks are in here
subset(links, peak %in% top_5_each$gene)
r_feature_plot(activity_ob, reduction = "umap_rna", features = "KCTD8") |
r_dim_plot(activity_ob, reduction = "umap_rna", group.by = c("RNA_cluster", "ATAC_cluster"))
table(activity_small$RNA_cluster, activity_small$ATAC_cluster)

##
# NO LINKED PEAKS FOUND IN OBJECT!! HOW?????
# ANSWER: YOUR DUMB ASS WASN'T USING THE peak COLUMN FROM THE LINKPEAKS OUTPUT
##

#lets try closest feature instead and filter against degs
close_genes <- ClosestFeature(activity_small, regions = diff_peaks$gene)
head(close_genes)

subset(close_genes,
       gene_name %in% degs$gene &
       distance == 0) %>%
    pull(gene_name)

#motif analysis
