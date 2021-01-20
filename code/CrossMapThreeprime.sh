#!/bin/bash

#SBATCH --job-name=CrossmapChimp3prime
#SBATCH --time=24:00:00
#SBATCH --partition=broadwl
#SBATCH --mem=36G
#SBATCH --tasks-per-node=4
#SBATCH --mail-type=END
#SBATCH --output=CrossmapChimp3prime.out
#SBATCH --error=CrossmapChimp3prime.err



source ~/activate_anaconda.sh
conda activate comp_threeprime_env



for i in $(ls ../Chimp/data/sort_clean/*.bam)
do
describer=$(echo ${i} | cut -d/ -f5| cut -d_ -f3-4 | cut -d- -f1)
echo $describer
CrossMap.py bam ../data/chainFiles/panTro6ToHg38.over.chain $i ../Chimp/data/sort_clean_hg38/$describer-hg38.bam -a
done
