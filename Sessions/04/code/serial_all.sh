#!/bin/bash
#SBATCH --mail-user=first.last@nationwidechildrens.org
#SBATCH --mail-type=FAIL,REQUEUE
#SBATCH --output=slurm-serial-all-%j.out
#SBATCH --cpus-per-task=1
#SBATCH --time=0-01:00:00
#SBATCH --mem-per-cpu=32G

set -e
ml purge
ml load STAR/2.7.9a
for((i=0; i<=4; i++)); do 
    STAR --outSAMtype BAM Unsorted --runThreadN 1 --outSAMstrandField intronMotif --genomeDir /reference/homo_sapiens/GRCh38/ensembl/release-96/Sequence/STARIndex/2.7.9/ --outFileNamePrefix parallel_all_${i} --readFilesIn ../benchmark_ex/data/sample${i}.fq
done
