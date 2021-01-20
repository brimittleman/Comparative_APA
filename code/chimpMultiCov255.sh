#!/bin/bash

#SBATCH --job-name=chimpMultiCov255
#SBATCH --time=24:00:00
#SBATCH --partition=broadwl
#SBATCH --mem=16G
#SBATCH --mail-type=END
#SBATCH --output=chimpMultiCov.out
#SBATCH --error=chimpMultiCov.err



source ~/activate_anaconda.sh
conda activate comp_threeprime_env
#opp strand

bedtools multicov -q 255 -S  -bams ../Chimp/data/sort_clean/*_N-clean.sort.bam -bed ../data/PAS_doubleFilter/PAS_doublefilter_either_ChimpCoordChimpUsage.sort.bed > ../data/testQuant/Chimp_DF_PAS_255.txt
