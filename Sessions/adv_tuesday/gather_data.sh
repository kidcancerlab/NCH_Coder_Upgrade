#!/usr/bin/bash

# Download the data to be used with the advanced Tuesday session on single cell ATAC-seq analysis

cd /home/gdworkshop/lab/Sessions/adv_tuesday

# https://www.10xgenomics.com/datasets/8k-adult-mouse-cortex-cells-atac-v2-chromium-controller-2-standard
# Single Cell ATAC Dataset by Cell Ranger ATAC 2.1.0
# Adult mouse cortex was obtained by 10x Genomics from BrainBits, LLC (PN C57ACX). Tissue was dissociated and nuclei were isolated using the demonstrated protocol Nuclei Isolation from Mouse Brain Tissue for Single Cell ATAC Sequencing (CG000212).

# ATAC libraries were generated as described in the Chromium Single Cell ATAC Reagent Kits User Guide (v2 chemistry) (CG000496) using the Chromium Controller and sequenced on Illumina NovaSeq 6000 to approximately 60k read pairs per cell. 7,729 cortex nuclei were recovered.

# Paired-end, dual indexing:

# 50 cycles Read 1
# 8 cycles i7 (sample index)
# 16 cycles i5 (10x barcode)
# 49 cycles Read 2

## Output Files
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller/8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller_analysis.tar.gz
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller/8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller_filtered_peak_bc_matrix.tar.gz
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller/8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller_filtered_peak_bc_matrix.h5
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller/8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller_filtered_tf_bc_matrix.tar.gz
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller/8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller_filtered_tf_bc_matrix.h5
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller/8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller_fragments.tsv.gz
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller/8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller_fragments.tsv.gz.tbi
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller/8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller_peak_annotation.tsv
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller/8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller_peak_motif_mapping.bed
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller/8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller_peaks.bed
wget https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-atac/2.1.0/8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller/8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller_possorted_bam.bam
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller/8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller_possorted_bam.bam.bai
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller/8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller_raw_peak_bc_matrix.tar.gz
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller/8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller_raw_peak_bc_matrix.h5
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller/8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller_singlecell.csv
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller/8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller_summary.csv
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller/8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller_summary.json
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller/8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller_web_summary.html
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller/8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller_cloupe.cloupe


# https://www.10xgenomics.com/datasets/8k-adult-mouse-cortex-cells-atac-v2-chromium-x-2-standard
# Single Cell ATAC Dataset by Cell Ranger ATAC 2.1.0
# Adult mouse cortex was obtained by 10x Genomics from BrainBits, LLC (PN C57ACX). Tissue was dissociated and nuclei were isolated using the demonstrated protocol Nuclei Isolation from Mouse Brain Tissue for Single Cell ATAC Sequencing (CG000212).

# ATAC libraries were generated as described in the Chromium Single Cell ATAC Reagent Kits User Guide (v2 chemistry) (CG000496) using the Chromium X and sequenced on Illumina NovaSeq 6000 to approximately 55k read pairs per cell. 8,067 cortex nuclei were recovered.

# Paired-end, dual indexing:

# 50 cycles Read 1
# 8 cycles i7 (sample index)
# 16 cycles i5 (10x barcode)
# 49 cycles Read 2

# Output Files
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/8k_mouse_cortex_ATACv2_nextgem_Chromium_X/8k_mouse_cortex_ATACv2_nextgem_Chromium_X_analysis.tar.gz
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/8k_mouse_cortex_ATACv2_nextgem_Chromium_X/8k_mouse_cortex_ATACv2_nextgem_Chromium_X_filtered_peak_bc_matrix.tar.gz
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/8k_mouse_cortex_ATACv2_nextgem_Chromium_X/8k_mouse_cortex_ATACv2_nextgem_Chromium_X_filtered_peak_bc_matrix.h5
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/8k_mouse_cortex_ATACv2_nextgem_Chromium_X/8k_mouse_cortex_ATACv2_nextgem_Chromium_X_filtered_tf_bc_matrix.tar.gz
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/8k_mouse_cortex_ATACv2_nextgem_Chromium_X/8k_mouse_cortex_ATACv2_nextgem_Chromium_X_filtered_tf_bc_matrix.h5
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/8k_mouse_cortex_ATACv2_nextgem_Chromium_X/8k_mouse_cortex_ATACv2_nextgem_Chromium_X_fragments.tsv.gz
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/8k_mouse_cortex_ATACv2_nextgem_Chromium_X/8k_mouse_cortex_ATACv2_nextgem_Chromium_X_fragments.tsv.gz.tbi
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/8k_mouse_cortex_ATACv2_nextgem_Chromium_X/8k_mouse_cortex_ATACv2_nextgem_Chromium_X_peak_annotation.tsv
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/8k_mouse_cortex_ATACv2_nextgem_Chromium_X/8k_mouse_cortex_ATACv2_nextgem_Chromium_X_peak_motif_mapping.bed
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/8k_mouse_cortex_ATACv2_nextgem_Chromium_X/8k_mouse_cortex_ATACv2_nextgem_Chromium_X_peaks.bed
wget https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-atac/2.1.0/8k_mouse_cortex_ATACv2_nextgem_Chromium_X/8k_mouse_cortex_ATACv2_nextgem_Chromium_X_possorted_bam.bam
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/8k_mouse_cortex_ATACv2_nextgem_Chromium_X/8k_mouse_cortex_ATACv2_nextgem_Chromium_X_possorted_bam.bam.bai
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/8k_mouse_cortex_ATACv2_nextgem_Chromium_X/8k_mouse_cortex_ATACv2_nextgem_Chromium_X_raw_peak_bc_matrix.tar.gz
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/8k_mouse_cortex_ATACv2_nextgem_Chromium_X/8k_mouse_cortex_ATACv2_nextgem_Chromium_X_raw_peak_bc_matrix.h5
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/8k_mouse_cortex_ATACv2_nextgem_Chromium_X/8k_mouse_cortex_ATACv2_nextgem_Chromium_X_singlecell.csv
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/8k_mouse_cortex_ATACv2_nextgem_Chromium_X/8k_mouse_cortex_ATACv2_nextgem_Chromium_X_summary.csv
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/8k_mouse_cortex_ATACv2_nextgem_Chromium_X/8k_mouse_cortex_ATACv2_nextgem_Chromium_X_summary.json
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/8k_mouse_cortex_ATACv2_nextgem_Chromium_X/8k_mouse_cortex_ATACv2_nextgem_Chromium_X_web_summary.html
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/8k_mouse_cortex_ATACv2_nextgem_Chromium_X/8k_mouse_cortex_ATACv2_nextgem_Chromium_X_cloupe.cloupe