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


Rscript ../../leafcutter/scripts/ds_plots.R  -e ../data/DiffSplice/hg38_ncbiRefseq_exonsfixed.gz ../data/DiffSplice_liftedJunc/MergeCombined_perind_numers.counts.fixed.gz ../data/DiffSplice_liftedJunc/groups_file.txt ../data/DiffSplice_liftedJunc/MergedRes_cluster_significance.txt -f 0.05 --output ../data/DiffSplice_liftedJunc/ds_plots.pdf
