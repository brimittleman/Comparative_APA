#!/bin/bash


#SBATCH --job-name=makeDIC
#SBATCH --output=makeDIC.out
#SBATCH --error=makeDIC.err
#SBATCH --time=10:00:00
#SBATCH --partition=broadwl
#SBATCH --mem=36G
#SBATCH --mail-type=END

source ~/activate_anaconda.sh
conda activate comp_threeprime_env

while read gene;
do echo "$gene"
Rscript PlotNuclearUsagebySpecies_DF_4DIC.R  -g "$gene"
done < ../data/DIC_Viz/DiCGenes_5fdr.txt
