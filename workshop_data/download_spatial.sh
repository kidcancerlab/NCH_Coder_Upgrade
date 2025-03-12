#!/bin/bash
my_loc=$(pwd)
set -e

if [ ! -d input/spatial ]
then
    mkdir -p input/spatial
fi

cd input/spatial/

# Download data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193107

# download GSE193107 data that has scRNA-seq, scATAC-seq and Visium data
wget \
    -O GSE193107_RAW.tar \
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE193107&format=file"

# decompress data
tar -xvf GSE193107_RAW.tar
rm GSE193107_RAW.tar

for s in `ls *.tar.gz`; do
    sample=`basename ${s} .tar.gz`
    echo ${sample}
    mkdir ${sample}
    mv ${sample}.tar.gz ${sample}.tif.gz ${sample}/
    cd ${sample}
    tar --strip-components=6 -zxvf ${sample}.tar.gz
    rm ${sample}.tar.gz
    mkdir spatial
    mv *.json *.png *.csv spatial
    cp spatial/tissue_hires_image.png spatial/tissue_lowres_image.png
    cd ../
done

cd ${my_loc}
