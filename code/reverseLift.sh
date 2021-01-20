#!/bin/bash

#SBATCH --job-name=revLift
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=revLift.out
#SBATCH --error=revLift.err
#SBATCH --partition=broadwl
#SBATCH --mem=36G
#SBATCH --mail-type=END

source ~/activate_anaconda.sh
conda activate comp_threeprime_env


#liftover  oldFile map.chain newFile unMapped

liftOver ../data/primaryLift/chimp_APApeaks_primarylift2Human.bed  ../data/chainFiles/hg38ToPanTro6.over.chain  ../data/reverseLift/chimp_APApeaks_primarylift2Human_rev2Chimp.bed ../data/reverseLift/chimp_APApeaks_primarylift2Human_rev2Chimp_UNLIFTED.bed


liftOver ../data/primaryLift/human_APApeaks_primarylift2Chimp.bed  ../data/chainFiles/panTro6ToHg38.over.chain ../data/reverseLift/human_APApeaks_primarylift2Chimp_rev2Human.bed ../data/reverseLift/human_APApeaks_primarylift2Human_rev2Human_UNLIFTED.bed
