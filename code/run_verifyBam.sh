#!/bin/bash


#SBATCH --job-name=run_verifybam
#SBATCH --output=run_verifybam.out
#SBATCH --error=run_verifybam.err
#SBATCH --time=36:00:00
#SBATCH --partition=broadwl
#SBATCH --constraint=edr
#SBATCH --nodes=6
#SBATCH --ntasks-per-node=13
#SBATCH --mem=36G
#SBATCH --mail-type=END

source ~/activate_anaconda.sh
conda activate comp_threeprime_env


#/project2/gilad/reem/vcf_fromjohn/sortedsnps.hg38liftover.exons.vcf.gz
#../Human/data/sort_clean/human_combined_18510_T-clean.sort.bam


i=$1

describer=$(echo ${i} | sed -e 's/.*\human_combined_//' | sed -e "s/clean.sort.bam$//")
verifyBamID --vcf /project2/gilad/reem/vcf_fromjohn/sortedsnps.hg38liftover.exons.vcf.gz --bam ${i} --best --ignoreRG --out ../Human/data/verifyBam/${describer}.verify
