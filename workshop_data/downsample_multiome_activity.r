library(Seurat)
library(DropletUtils) #enables us to write10xCounts
library(tidyverse)
library(Signac)
old_data <- Read10X_h5("filtered_feature_bc_matrix.h5")

#get barcodes to keep
set.seed(17)
bcs_keep <- sample(colnames(old_data[["Gene Expression"]]), size = 2000)

new_gene_exp <- old_data[["Gene Expression"]][, bcs_keep]

#write off h5 with atac data too
new_atac <- old_data[["Peaks"]][, bcs_keep]

new_data <- list("Gene Expression" = new_gene_exp,
                 "Peaks" = new_atac)

write.table(bcs_keep,
            "../../downsampled_data/downsampled_bcs.tsv",
            row.names = FALSE,
            quote = FALSE,
            col.names = FALSE)
qs::qsave(new_data, "../../downsampled_data/gene_exp_peaks_list.qs")