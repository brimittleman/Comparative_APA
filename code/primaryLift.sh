#!/bin/bash

#SBATCH --job-name=primaryLift
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=primaryLift.out
#SBATCH --error=primaryLift.err
#SBATCH --partition=broadwl
#SBATCH --mem=36G
#SBATCH --mail-type=END

source ~/activate_anaconda.sh
conda activate comp_threeprime_env


#liftover  oldFile map.chain newFile unMapped

liftOver ../data/cleanPeaks_byspecies/human_APApeaks.ALLChrom.Filtered.Named.Cleaned_100bpreg.bed   ../data/chainFiles/hg38ToPanTro6.over.chain  ../data/primaryLift/human_APApeaks_primarylift2Chimp.bed  ../data/primaryLift/human_APApeaks_primarylift2Chimp_UNLIFTED.bed


liftOver ../data/cleanPeaks_byspecies/chimp_APApeaks.ALLChrom.Filtered.Named.Cleaned_100bpreg.bed ../data/chainFiles/panTro6ToHg38.over.chain ../data/primaryLift/chimp_APApeaks_primarylift2Human.bed  ../data/primaryLift/chimp_APApeaks_primarylift2Human_UNLIFTED.bed
