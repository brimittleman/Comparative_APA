#!/bin/bash

#SBATCH --job-name=Lift5perPASbed
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=Lift5perPASbed.out
#SBATCH --error=Lift5perPASbed.err
#SBATCH --partition=broadwl
#SBATCH --mem=36G
#SBATCH --mail-type=END

source ~/activate_anaconda.sh
conda activate comp_threeprime_env


# oldFile map.chain newFile unMapped


liftOver ../data/PAS/PAS_5perc_either_HumanCoordChimpUsage.bed  ../data/chainFiles/hg38ToPanTro6.over.chain ../data/PAS/PAS_5perc_either_ChimpCoordChimpUsage.bed  ../data/PAS/PAS_5perc_either_unlifted
