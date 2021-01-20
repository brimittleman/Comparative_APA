#!/bin/bash

#SBATCH --job-name=recChimpback2Human
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=recChimpback2Human.out
#SBATCH --error=recChimpback2Human.err
#SBATCH --partition=broadwl
#SBATCH --mem=36G
#SBATCH --mail-type=END

source ~/activate_anaconda.sh
conda activate comp_threeprime_env



liftOver  ../data/cleanPeaks_lifted/Chimp_PASregions.bed ../data/chainFiles/panTro6ToHg38.over.chain  ../data/cleanPeaks_lifted/Chimp_PASregions_humanCoord.bed ../data/cleanPeaks_lifted/Chimp_PASregions_humanCoord_notlift.bed
