#!/bin/bash
#SBATCH --mail-user=first.last@nationwidechildrens.org
#SBATCH --mail-type=FAIL,REQUEUE
#SBATCH --output=slurm-serial-single-5cpu-%j.out
#SBATCH --cpus-per-task=5
#SBATCH --time=0-01:00:00

set -e
ml purge
ml load STAR/2.7.9a
time STAR --outSAMtype BAM Unsorted --runThreadN 5 --outSAMstrandField intronMotif --genomeDir /reference/homo_sapiens/GRCh38/ensembl/release-96/Sequence/STARIndex/2.7.9/ --outFileNamePrefix serial0_cpu5 --readFilesIn ../benchmark_ex/data/sample0.fq
