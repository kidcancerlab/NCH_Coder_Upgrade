#!/bin/bash
my_loc=$(pwd)

set -e

if [ ! -d input/scRNA ]
then
    mkdir -p input/scRNA
fi

# We will store the raw data within input/scRNA
cd input/scRNA/

# Download data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129788
# paper: https://www.nature.com/articles/s41593-019-0491-3

# download GSE222510 data that has scRNA-seq
wget \
    -O GSE222510_RAW.tar \
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE222510&format=file"

# decompress data
tar -xvf GSE222510_RAW.tar
rm GSE222510_RAW.tar

# We don't need most of these samples. Keeping only old and young
rm GSM*OY*
rm GSM*YO*
rm GSM*YY*
rm GSM*OO*

kept_samples=(GSM6925133_OX1X
              GSM6925134_OX2X
              GSM6925135_OX3X
              GSM6925136_OX4X
              GSM6925137_OX5X
              GSM6925138_OX6X
              GSM6925139_OX7X
              GSM6925140_OX8X
              GSM6925161_YX1L
              GSM6925162_YX2L
              GSM6925163_YX3R
              GSM6925164_YX4R
              GSM6925165_YX5R
              GSM6925166_YX6L
              GSM6925167_YX7R
              GSM6925168_YX8L)

for this_sample in ${kept_samples[@]}
do
    mkdir ${this_sample}
    mv ${this_sample}_* ${this_sample}/
    rename ${this_sample}_ "" ${this_sample}/*

    # The file is improperly formatted for Seurat to take
    zcat ${this_sample}/genes.tsv.gz \
        | perl -ne 'chomp; print $_, "\t", $_, "\tGene Expression\n"' \
        | gzip \
        > ${this_sample}/features.tsv

    rm ${this_sample}/genes.tsv.gz
done

wget \
    -O GSE222510_filtered_metadata.txt.gz \
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE222nnn/GSE222510/suppl/GSE222510%5Ffiltered%5Fmetadata.txt.gz"

cd ${my_loc}
