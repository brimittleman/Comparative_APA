#!/bin/bash

#SBATCH --job-name=mergeandsort_ChimpinHuman
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=mergeandsort_ChimpinHuman.out
#SBATCH --error=mergeandsort_ChimpinHuman.err
#SBATCH --partition=bigmem2
#SBATCH --mem=100G
#SBATCH --mail-type=END

source ~/activate_anaconda.sh
conda activate comp_threeprime_env


samtools merge ../Chimp/data/RNAseq/mergedBAM_hg38/Chimp_RNAMerged_inhuman.bam  ../Chimp/data/RNAseq/Sort_hg38/*.bam.sorted.bam

samtools sort ../Chimp/data/RNAseq/mergedBAM_hg38/Chimp_RNAMerged_inhuman.bam -o ../Chimp/data/RNAseq/mergedBAM_hg38/Chimp_RNAMerged_inhuman.sort.bam 

samtools index  ../Chimp/data/RNAseq/mergedBAM_hg38/Chimp_RNAMerged_inhuman.sort.bam
