#!/bin/bash

#SBATCH --job-name=buildLeafviz
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=buildLeafviz.out
#SBATCH --error=buildLeafviz.err
#SBATCH --partition=broadwl
#SBATCH --mem=36G
#SBATCH --mail-type=END

module load Anaconda3/5.3.0
source  activate leafcutter


/project2/gilad/briana/leafcutter/leafviz/prepare_results.R -o ../data/leafviz/ChimpvHumanLeafviz -m ../data/DiffSplice_liftedJunc/groups_file.txt -c LeafvizHC  ../data/DiffSplice_liftedJunc/MergeCombined_perind_numers.counts.fixed.gz  ../data/DiffSplice_liftedJunc/MergedRes_cluster_significance.txt ../data/DiffSplice_liftedJunc/MergedRes_effect_sizes.txt ../data/leafviz/hg38 
