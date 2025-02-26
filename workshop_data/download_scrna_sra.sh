#!/bin/bash

set -e

if [ ! -d tempdir_fastqs ]
then
    mkdir -p tempdir_fastqs/fastqs
fi

sra_numbers=(SRR8895023
             SRR8895024
             SRR8895025
             SRR8895026
             SRR8895027
             SRR8895028
             SRR8895029
             SRR8895030
             SRR8895031
             SRR8895032
             SRR8895033
             SRR8895034
             SRR8895035
             SRR8895036
             SRR8895037
             SRR8895038)

cd tempdir_fastqs

parallel
    -j 10 \
    "prefetch \
        -C yes \
        --max-size 300g \
        {}" \
    ::: ${sra_numbers[@]}

# Make a directory for each SRA
for this_sra in ${sra_numbers[@]}
do
    mkdir fastqs/${this_sra}
done

parallel \
    -j 10 \
    "fasterq-dump \
        --split-files \
        --outdir fastqs/{} \
        {}" \
    ::: ${sra_numbers[@]}

# gzip the fastqs
parallel \
    "pigz {}" \
    ::: fastqs/${this_sra}/*.fastq

# Need to make the fastqs match the illumina file naming convention
# sample_name_S1_L00[1-0]_R[12]_001.fastq.gz
for this_sra in ${sra_numbers[@]}
do
    mv fastqs/${this_sra}/${this_sra}_1.fastq.gz fastqs/${this_sra}/${this_sra}_S1_L001_R1_001.fastq.gz
    mv fastqs/${this_sra}/${this_sra}_2.fastq.gz fastqs/${this_sra}/${this_sra}_S1_L001_R2_001.fastq.gz
done

# cleanup the sra files
for this_sra in ${sra_numbers[@]}
do
    rm -r ${this_sra}
done
