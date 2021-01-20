#!/bin/bash


#SBATCH --job-name=verifybam4973HumanMap
#SBATCH --output=verifybam4973HumanMap.out
#SBATCH --error=verifybam4973HumanMap.err
#SBATCH --time=36:00:00
#SBATCH --partition=broadwl
#SBATCH --ntasks-per-node=13
#SBATCH --mem=16G
#SBATCH --mail-type=END

source ~/activate_anaconda.sh
conda activate comp_threeprime_env




verifyBamID --vcf /project2/gilad/reem/vcf_fromjohn/sortedsnps.hg38liftover.exons.vcf.gz --bam ../data/TwoBadSampleAnalysis/NA4973inHuman-sort.bam --ignoreRG --best --out ../data/TwoBadSampleAnalysis/4973_inHuman.verify
