#!/bin/bash

#SBATCH --job-name=buildLeafviz_leafanno
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=buildLeafviz_leafanno.out
#SBATCH --error=buildLeafviz_leafanno.err
#SBATCH --partition=broadwl
#SBATCH --mem=36G
#SBATCH --mail-type=END

module load Anaconda3/5.3.0
source  activate leafcutter


/project2/gilad/briana/leafcutter/leafviz/prepare_results.R -o ../data/leafviz/ChimpvHumanLeafviz_leafanno -m ../data/DiffSplice_liftedJunc/groups_file.txt -c LeafvizHC_leafanno  ../data/DiffSplice_liftedJunc/MergeCombined_perind_numers.counts.fixed.gz  ../data/DiffSplice_liftedJunc/MergedRes_cluster_significance.txt ../data/DiffSplice_liftedJunc/MergedRes_effect_sizes.txt ../data/leafviz/annotation_codes/gencode_hg38/gencode_hg38
