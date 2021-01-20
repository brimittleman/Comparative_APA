#!/bin/bash

#SBATCH --job-name=humanMultiCovInclusive
#SBATCH --time=24:00:00
#SBATCH --partition=broadwl
#SBATCH --mem=32G
#SBATCH --mail-type=END
#SBATCH --output=humanMultiCovInclusive.out
#SBATCH --error=humanMultiCovInclusive.err



source ~/activate_anaconda.sh
conda activate comp_threeprime_env
#actual strand like in original

bedtools multicov -s -split -bams ../Human/data/sort_clean/*_N-clean.sort.bam -bed ../data/testQuant/Human_ALLpeaks.Bothstrand.bed > ../data/testQuant/Human_AllPeaksInclusive.txt
