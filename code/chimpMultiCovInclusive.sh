#!/bin/bash

#SBATCH --job-name=chimpMultiCovInclusive
#SBATCH --time=24:00:00
#SBATCH --partition=broadwl
#SBATCH --mem=32G
#SBATCH --mail-type=END
#SBATCH --output=chimpMultiCovInclusive.out
#SBATCH --error=chimpMultiCovInclusive.err



source ~/activate_anaconda.sh
conda activate comp_threeprime_env
#actual strand like in original

bedtools multicov -s -split  -bams ../Chimp/data/sort_clean/*_N-clean.sort.bam -bed ../data/testQuant/Chimp_ALLpeaks.Bothstrand.bed > ../data/testQuant/Chimp_AllPeaksInclusive.txt
