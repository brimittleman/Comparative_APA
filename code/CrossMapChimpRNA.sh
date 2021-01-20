#!/bin/bash

#SBATCH --job-name=CrossmapChimpRNA
#SBATCH --time=24:00:00
#SBATCH --partition=broadwl
#SBATCH --mem=36G
#SBATCH --tasks-per-node=4
#SBATCH --mail-type=END
#SBATCH --output=CrossmapChimpRNA.out
#SBATCH --error=CrossmapChimpRNA.err



source ~/activate_anaconda.sh
conda activate comp_threeprime_env



for i in $(ls ../Chimp/data/RNAseq/sort/*.bam)
do
describer=$(echo ${i} | cut -d/ -f6 | cut -d- -f2-5)
echo $describer
CrossMap.py bam ../data/chainFiles/panTro6ToHg38.over.chain $i ../Chimp/data/RNAseq/Sort_hg38/$describer-hg38.bam -a
done
