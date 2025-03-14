#!/bin/bash
my_loc=$(pwd)

set -e

if [ ! -d input/multiomics ]
then
    mkdir -p input/multiomics
fi

cd input/multiomics

# From 10x genomics website
# https://www.10xgenomics.com/datasets/frozen-human-healthy-brain-tissue-3-k-1-standard-2-0-0

# Output Files
# wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_cloupe.cloupe
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_analysis.tar.gz
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_web_summary.html
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_summary.csv
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_per_barcode_metrics.csv
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_filtered_feature_bc_matrix.tar.gz
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_filtered_feature_bc_matrix.h5
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_raw_feature_bc_matrix.tar.gz
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_raw_feature_bc_matrix.h5
# wget https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_gex_possorted_bam.bam
# wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_gex_possorted_bam.bam.bai
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_gex_molecule_info.h5
# wget https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_atac_possorted_bam.bam
# wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_atac_possorted_bam.bam.bai
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_atac_fragments.tsv.gz
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_atac_fragments.tsv.gz.tbi
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_atac_peaks.bed
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_atac_cut_sites.bigwig
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_atac_peak_annotation.tsv

cd ${my_loc}
