#!/bin/bash

#SBATCH --job-name=hg19MapStats
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=hg19MapStats.out
#SBATCH --error=hg19MapStats.err
#SBATCH --partition=broadwl
#SBATCH --mem=36G
#SBATCH --mail-type=END

source ~/activate_anaconda.sh
conda activate comp_threeprime_env




touch "/project2/gilad/briana/Comparative_APA/data/MapStats/Human_hg19_RNAseq_mapstats.txt"


for i in $(ls /project2/gilad/briana/Comparative_APA/Human/data/RNAseq/bam_hg19/*new-sort.bam)
do
echo $i
samtools view $i  | cut -f 2,5,12 >> /project2/gilad/briana/Comparative_APA/data/MapStats/Human_hg19_RNAseq_mapstats.txt
done
