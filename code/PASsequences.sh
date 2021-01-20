#!/bin/bash

#SBATCH --job-name=PAS_sequence
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=PAS_sequence.out
#SBATCH --error=PAS_sequence.err
#SBATCH --partition=broadwl
#SBATCH --mem=36G
#SBATCH --mail-type=END

source ~/activate_anaconda.sh
conda activate comp_threeprime_env


bedtools nuc -s -seq  -pattern "AATAAA" -C -fi /project2/gilad/kenneth/References/human/genome/hg38.fa -bed ../data/PAS/PAS_5perc_either_HumanCoordHummanUsage.sort.bed > ../data/SignalSites/PAS_5perc_either_HumanCoordHummanUsage_nuc.txt

bedtools nuc -s -seq  -pattern "AATAAA" -C -fi /project2/gilad/briana/genome_anotation_data/Chimp_genome/panTro6.fa -bed ../data/PAS/PAS_5perc_either_ChimpCoordChimpUsage.sort.bed > ../data/SignalSites/PAS_5perc_either_ChimpCoordChimpUsage_nuc.txt
