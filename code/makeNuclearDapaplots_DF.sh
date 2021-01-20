#!/bin/bash


#SBATCH --job-name=makeNuclearPlotsDF
#SBATCH --output=makeNuclearPlotsDF.out
#SBATCH --error=makeNuclearPlotsDF.err
#SBATCH --time=10:00:00
#SBATCH --partition=broadwl
#SBATCH --mem=36G
#SBATCH --mail-type=END

source ~/activate_anaconda.sh
conda activate comp_threeprime_env

while read gene;
do echo "$gene"
Rscript PlotNuclearUsagebySpecies_DF.R -g "$gene"
done < ../data/DiffIso_Nuclear_DF/SignifianceEitherGENES_Nuclear_noHead.txt
