#!/bin/bash

#SBATCH --job-name=ChimpchromOrder
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=chromOrder.out
#SBATCH --error=chromOrder.err
#SBATCH --partition=broadwl
#SBATCH --mem=36G
#SBATCH --mail-type=END


cat /project2/gilad/briana/Comparative_APA/Chimp/data/bed_clean/* | cut -f1 | sort | uniq > /project2/gilad/briana/genome_anotation_data/Pantro5/panTro5.chromOrder2
