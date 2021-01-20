#!/bin/bash

#SBATCH --job-name=overlapPAS
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=overlapPAS.out
#SBATCH --error=overlapPAS.err
#SBATCH --partition=broadwl
#SBATCH --mem=36G
#SBATCH --mail-type=END

source ~/activate_anaconda.sh
conda activate comp_threeprime_env

bedtools intersect -f 0.1 -r -s -wo -sorted -a ../data/cleanPeaks_lifted/Human_PASregions.sort.bed -b ../data/cleanPeaks_lifted/Chimp_PASregions_humanCoord.sort.bed > ../data/cleanPeaks_lifted/PASregions_identifiedbothTEST.txt

bedtools intersect -f 0.625 -r -s -wo -sorted -a ../data/cleanPeaks_lifted/Human_PASregions.sort.bed -b ../data/cleanPeaks_lifted/Chimp_PASregions_humanCoord.sort.bed > ../data/cleanPeaks_lifted/PASregions_identifiedboth.txt


#only human
bedtools intersect -f 0.625 -r -s -sorted -v -a ../data/cleanPeaks_lifted/Human_PASregions.sort.bed -b ../data/cleanPeaks_lifted/Chimp_PASregions_humanCoord.sort.bed > ../data/cleanPeaks_lifted/PASregions_identifiedHuman.txt



#only chimp
bedtools intersect -f 0.625 -r -s -sorted -v -a ../data/cleanPeaks_lifted/Chimp_PASregions_humanCoord.sort.bed -b ../data/cleanPeaks_lifted/Human_PASregions.sort.bed > ../data/cleanPeaks_lifted/PASregions_identifiedChimp.txt
