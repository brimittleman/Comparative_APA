#!/bin/bash

#SBATCH --job-name=humanMultiCov
#SBATCH --time=24:00:00
#SBATCH --partition=broadwl
#SBATCH --mem=16G
#SBATCH --mail-type=END
#SBATCH --output=humanMultiCov.out
#SBATCH --error=humanMultiCov.err



source ~/activate_anaconda.sh
conda activate comp_threeprime_env
#opp strand 

bedtools multicov -S  -bams ../Human/data/sort_clean/*_N-clean.sort.bam -bed ../data/PAS_doubleFilter/PAS_doublefilter_either_HumanCoordHummanUsage.sort.bed > ../data/testQuant/Human_DF_PAS.txt
