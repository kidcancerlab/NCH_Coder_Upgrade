#!/bin/bash

set -e

my_loc=`pwd`

if [ ! -d input/scRNA ]
then
    mkdir -p input/scRNA
fi

# We will store the raw data within input/scRNA
cd input/scRNA/

# Download data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129788

# download GSE129788 data that has scRNA-seq
wget \
    -O GSE129788_RAW.tar \
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE129788&format=file"

# decompress data
tar -xvf GSE129788_RAW.tar
rm GSE129788_RAW.tar

wget \
    -O GSE129788_Supplementary_meta_data_Cell_Types_Etc.txt.gz \
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE129788&format=file&file=GSE129788%5FSupplementary%5Fmeta%5Fdata%5FCell%5FTypes%5FEtc%2Etxt%2Egz"

cd ${my_loc}
