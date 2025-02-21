#!/bin/bash

home="/igm/home/cnh008/projects/scrgot/scrgot_2025/mouse_brain/spatial/raw"
cd ${home}

# Download data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193107

# download GSE193107 data that has scRNA-seq, scATAC-seq and Visium data
wget --recursive --no-parent -nd https://ftp.ncbi.nlm.nih.gov/geo/series/GSE193nnn/GSE193107/suppl/

# decompress data
tar -xvf GSE193107_RAW.tar

for s in `ls *.tar.gz`; do
    sample=`basename ${s} .tar.gz`
    mkdir ${sample}
    mv ${sample}.tar.gz ${sample}.tif.gz ${sample}
    cd ${sample}
    tar --strip-components=6 -zxvf ${sample}.tar.gz
    mkdir spatial
    mv *.json *.png *.csv spatial
    cp spatial/tissue_hires_image.png spatial/tissue_lowres_image.png
    cd ${home}
done