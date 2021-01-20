#!/bin/bash


#SBATCH --job-name=verifybam4973
#SBATCH --output=verifybam4973.out
#SBATCH --error=verifybam4973.err
#SBATCH --time=36:00:00
#SBATCH --partition=broadwl
#SBATCH --ntasks-per-node=13
#SBATCH --mem=16G
#SBATCH --mail-type=END


source ~/activate_anaconda.sh
conda activate comp_threeprime_env

verifyBamID --vcf /project2/gilad/reem/vcf_fromjohn/sortedsnps.hg38liftover.exons.vcf.gz --bam /project2/gilad/briana/Comparative_APA/Chimp/data/sort_clean/chimp_combined_4973_N-clean.sort.bam --ignoreRG --best --out ../Human/data/verifyBam/4973_N.verify

verifyBamID --vcf /project2/gilad/reem/vcf_fromjohn/sortedsnps.hg38liftover.exons.vcf.gz --bam /project2/gilad/briana/Comparative_APA/Chimp/data/sort_clean/chimp_combined_4973_T-clean.sort.bam --ignoreRG --best --out ../Human/data/verifyBam/4973_T.verify
