#!/bin/bash

#SBATCH --job-name=IntersectPASandOrtho
#SBATCH --time=24:00:00
#SBATCH --partition=broadwl
#SBATCH --mem=16G
#SBATCH --mail-type=END
#SBATCH --output=IntersectPASandOrtho.out
#SBATCH --error=IntersectPASandOrtho.err



source ~/activate_anaconda.sh
conda activate comp_threeprime_env


bedtools intersect -sorted -s  -a ../data/PAS_doubleFilter/PAS_doublefilter_either_HumanCoordHummanUsage.bed -b ../data/OrthoExonBed/human.noM.sort.bed > ../data/OrthoExonBed/allPASinOrtho.bed

bedtools intersect -sorted -s -v -a ../data/PAS_doubleFilter/PAS_doublefilter_either_HumanCoordHummanUsage.bed -b ../data/OrthoExonBed/human.noM.sort.bed > ../data/OrthoExonBed/allPAS_NOT_inOrtho.bed
