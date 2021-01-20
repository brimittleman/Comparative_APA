#!/bin/bash

#SBATCH --job-name=NuclearPlotsDEandDiffDom_4
#SBATCH --output=NuclearPlotsDEandDiffDom_4.out
#SBATCH --error=NuclearPlotsDEandDiffDom_4.err
#SBATCH --time=10:00:00
#SBATCH --partition=broadwl
#SBATCH --mem=36G
#SBATCH --mail-type=END

source ~/activate_anaconda.sh
conda activate comp_threeprime_env

while read gene;
do echo "$gene"
Rscript PlotNuclearUsagebySpecies_DF_DEout.R  -g "$gene"
done < ../data/DiffDomandDE_example/genesfor4examples.txt
