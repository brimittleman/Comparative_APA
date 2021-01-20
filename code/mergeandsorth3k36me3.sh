#!/bin/bash

#SBATCH --job-name=mergeandsort_h3k36me3
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=mergeandsort_h3k36me3
#SBATCH --error=mergeandsort_h3k36me3
#SBATCH --partition=bigmem2
#SBATCH --mem=100G
#SBATCH --mail-type=END

source ~/activate_anaconda.sh
conda activate comp_threeprime_env

samtools merge ../data/H3K36me3/MergedIndiv.H3k36me3.bam  ../data/H3K36me3/*.bam

samtools sort ../data/H3K36me3/MergedIndiv.H3k36me3.bam  -o ../data/H3K36me3/MergedIndiv.H3k36me3.sort.bam 

samtools index ../data/H3K36me3/MergedIndiv.H3k36me3.sort.bam 
