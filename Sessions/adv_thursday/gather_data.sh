#!/usr/bin/bash

# Download the data to be used with the advanced Thursday session on single cell ATAC-seq analysis

cd /home/gdworkshop/lab/Sessions/adv_thursday

# https://www.10xgenomics.com/datasets/frozen-human-healthy-brain-tissue-3-k-1-standard-2-0-0
# Single Cell Multiome ATAC + Gene Expression Dataset by Cell Ranger ARC 2.0.0
# Flash-frozen human healthy brain tissue (cerebellum) was acquired from BioIVT Asterand®.
# Nuclei were isolated from a 2605mg section of flash-frozen healthy human brain following the Nuclei Isolation from Complex Tissues for Single Cell Multiome ATAC + Gene Expression Sequencing demonstrated protocol (CG000375). The isolated nuclei were stained with 7-AAD and then flow sorted using the BD FACSMelody™ cell sorter to clean up the nuclei from the debris. The 7-AAD positive sorted nuclei were then permeabilized using lysis buffer containing digitonin as per the demonstrated protocol.

# Sequencer: Novaseq 6000 v1 Kit (Forward Strand Dual-Index Workflow)
# Cycle numbers:
#   Gene Expression library: 28,10,10,90
#   ATAC library: 50,8,16,49

# Output Files
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_cloupe.cloupe
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_analysis.tar.gz
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_web_summary.html
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_summary.csv
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_per_barcode_metrics.csv
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_filtered_feature_bc_matrix.tar.gz
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_filtered_feature_bc_matrix.h5
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_raw_feature_bc_matrix.tar.gz
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_raw_feature_bc_matrix.h5
wget https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_gex_possorted_bam.bam
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_gex_possorted_bam.bam.bai
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_gex_molecule_info.h5
wget https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_atac_possorted_bam.bam
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_atac_possorted_bam.bam.bai
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_atac_fragments.tsv.gz
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_atac_fragments.tsv.gz.tbi
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_atac_peaks.bed
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_atac_cut_sites.bigwig
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_atac_peak_annotation.tsv