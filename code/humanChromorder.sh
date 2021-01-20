#!/bin/bash

#SBATCH --job-name=HumanchromOrder
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=HchromOrder.out
#SBATCH --error=HchromOrder.err
#SBATCH --partition=broadwl
#SBATCH --mem=36G
#SBATCH --mail-type=END


cat /project2/gilad/briana/Comparative_APA/Human/data/bed_clean/* | cut -f1 | sort | uniq > /project2/gilad/briana/genome_anotation_data/chromOrder.num.chrHuman2.txt
