#!/bin/bash


#SBATCH --job-name=run_Chimpverifybam
#SBATCH --output=run_Chimpverifybam.out
#SBATCH --error=run_Chimpverifybam.err
#SBATCH --time=36:00:00
#SBATCH --partition=broadwl
#SBATCH --constraint=edr
#SBATCH --nodes=6
#SBATCH --ntasks-per-node=13
#SBATCH --mem=36G
#SBATCH --mail-type=END

source ~/activate_anaconda.sh
conda activate comp_threeprime_env



i=$1

describer=$(echo ${i} | sed -e 's/.*chimp_combined_//'  | sed -e "s/-clean.sort.bam$//")
verifyBamID --vcf /project2/gilad/HumanChimpData/VCF/chimp.vcf.gz --bam ${i} --best --ignoreRG --out ../Chimp/data/verifyBam/${describer}.verify
