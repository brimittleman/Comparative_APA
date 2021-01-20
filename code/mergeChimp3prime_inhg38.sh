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

#total fraction

samtools merge ../Chimp/data/mergedbyFracBam_hg38/Chimp_TotalMerged_inhuman.bam  ../Chimp/data/sort_clean_hg38/*T*.sorted.bam

samtools sort ../Chimp/data/mergedbyFracBam_hg38/Chimp_TotalMerged_inhuman.bam -o ../Chimp/data/mergedbyFracBam_hg38/Chimp_TotalMerged_inhuman_sort.bam

samtools index ../Chimp/data/mergedbyFracBam_hg38/Chimp_TotalMerged_inhuman_sort.bam

#nuclear fractions

samtools merge ../Chimp/data/mergedbyFracBam_hg38/Chimp_NuclearMerged_inhuman.bam  ../Chimp/data/sort_clean_hg38/*N*.sorted.bam

samtools sort ../Chimp/data/mergedbyFracBam_hg38/Chimp_NuclearMerged_inhuman.bam -o  ../Chimp/data/mergedbyFracBam_hg38/Chimp_NuclearMerged_inhuman_sort.bam

samtools index ../Chimp/data/mergedbyFracBam_hg38/Chimp_NuclearMerged_inhuman_sort.bam
