#!/bin/bash

#SBATCH --job-name=MapStats
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=MapStats.out
#SBATCH --error=MapStats.err
#SBATCH --partition=broadwl
#SBATCH --mem=36G
#SBATCH --mail-type=END

source ~/activate_anaconda.sh
conda activate comp_threeprime_env




touch "/project2/gilad/briana/Comparative_APA/data/MapStats/Human_RNAseq_mapstats.txt"
touch "/project2/gilad/briana/Comparative_APA/data/MapStats/Chimp_RNAseq_mapstats.txt"

for i in $(ls /project2/gilad/briana/Comparative_APA/Human/data/RNAseq/sort/*.bam)
do
echo $i
samtools view $i  | cut -f 2,5,12 >> /project2/gilad/briana/Comparative_APA/data/MapStats/Human_RNAseq_mapstats.txt
done

for i in $(ls /project2/gilad/briana/Comparative_APA/Chimp/data/RNAseq/sort/*.bam)
do
echo $i
samtools view $i  | cut -f 2,5,12 >> /project2/gilad/briana/Comparative_APA/data/MapStats/Chimp_RNAseq_mapstats.txt
done
