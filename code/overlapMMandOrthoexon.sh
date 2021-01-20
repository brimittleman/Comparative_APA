#!/bin/bash

#SBATCH --job-name=IntersectMMandOrtho
#SBATCH --time=24:00:00
#SBATCH --partition=broadwl
#SBATCH --mem=16G
#SBATCH --mail-type=END
#SBATCH --output=IntersectMMandOrtho.out
#SBATCH --error=IntersectMMandOrtho.err



source ~/activate_anaconda.sh
conda activate comp_threeprime_env


bedtools intersect -sorted -s  -a ../data/multimap/allMM.sort.bed -b ../data/OrthoExonBed/human.noM.sort.bed > ../data/OrthoExonBed/allMMinOrtho.bed

bedtools intersect -sorted -s -v -a ../data/multimap/allMM.sort.bed -b ../data/OrthoExonBed/human.noM.sort.bed > ../data/OrthoExonBed/allMM_NOT_inOrtho.bed
