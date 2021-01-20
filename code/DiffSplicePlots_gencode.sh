#!/bin/bash

#SBATCH --job-name=DiffSplicePlots
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=DiffSplicePlots.out
#SBATCH --error=DiffSplicePlots.err
#SBATCH --partition=broadwl
#SBATCH --mem=36G
#SBATCH --mail-type=END

module load Anaconda3/5.3.0
source  activate leafcutter


Rscript ../../leafcutter/scripts/ds_plots.R  -e ../data/DiffSplice/gencode19_exons.txt.gz ../data/DiffSplice/BothSpec_perind.counts.gz ../data/DiffSplice/groups_file.txt ../data/DiffSplice/Gencode__cluster_significance.txt -f 0.05 --output ../data/DiffSplice/ds_plots_gencode.pdf
