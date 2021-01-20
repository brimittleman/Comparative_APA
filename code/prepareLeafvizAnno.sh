#!/bin/bash

#SBATCH --job-name=prepareAnnoLeafviz
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=prepareAnnoLeafviz.out
#SBATCH --error=prepareAnnoLeafviz.err
#SBATCH --partition=broadwl
#SBATCH --mem=36G
#SBATCH --mail-type=END

module load Anaconda3/5.3.0
source  activate leafcutter


/project2/gilad/briana/leafcutter/leafviz/gtf2leafcutter.pl -o ../data/leafviz/hg38  /project2/gilad/briana/genome_anotation_data/hg38_try2/hg38.gtf
