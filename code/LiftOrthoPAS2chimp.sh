#!/bin/bash

#SBATCH --job-name=LiftorthoPAS
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=LiftorthoPASt.out
#SBATCH --error=LiftorthoPAS.err
#SBATCH --partition=broadwl
#SBATCH --mem=36G
#SBATCH --mail-type=END

source ~/activate_anaconda.sh
conda activate comp_threeprime_env


#liftover  oldFile map.chain newFile unMapped

liftOver ../data/cleanPeaks_anno/AllPAS_postLift.sort_LocAnnoPARSED.bed    ../data/chainFiles/hg38ToPanTro6.over.chain  ../data/cleanPeaks_anno/AllPAS_postLift.sort_LocAnnoPARSED_chimpLoc.bed  ../data/cleanPeaks_anno/AllPAS_postLift.sort_LocAnnoPARSED_chimpLocUNLIFTED.bed
