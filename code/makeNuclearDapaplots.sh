#!/bin/bash


#SBATCH --job-name=makeNuclearPlots
#SBATCH --output=makeNuclearPlots.out
#SBATCH --error=makeNuclearPlots.err
#SBATCH --time=10:00:00
#SBATCH --partition=broadwl
#SBATCH --mem=36G
#SBATCH --mail-type=END

source ~/activate_anaconda.sh
conda activate comp_threeprime_env

while read gene;
do echo "$gene"
Rscript PlotNuclearUsagebySpecies.R -g "$gene"
done < ../data/DiffIso_Nuclear/SignifianceEitherGENES_Nuclear_noHead.txt
