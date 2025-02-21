#!/bin/bash

home="/igm/home/cnh008/projects/scrgot/scrgot_2025/mouse_brain/rna/raw"
cd ${home}

# Download data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129788

# download GSE193107 data that has scRNA-seq, scATAC-seq and Visium data
wget --recursive --no-parent -nd https://ftp.ncbi.nlm.nih.gov/geo/series/GSE129nnn/GSE129788/suppl/

# decompress data
tar -xvf GSE129788_RAW.tar