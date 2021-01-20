#!/bin/bash

#SBATCH --job-name=classifyLeafviz
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=classifyLeafviz.out
#SBATCH --error=classifyLeafviz.err
#SBATCH --partition=broadwl
#SBATCH --mem=36G
#SBATCH --mail-type=END

module load Anaconda3/5.3.0
source  activate leafcutter


Rscript ../data/leafviz/classify_clusters.R -o ../data/leafviz/ClassifiedResults ../data/leafviz/ChimpvHumanLeafviz_leafanno


Rscript ../data/leafviz/classify_clusters.R -o ../data/leafviz/ClassifiedResults_myAnno ../data/leafviz/ChimpvHumanLeafviz
