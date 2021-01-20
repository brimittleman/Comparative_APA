#!/bin/bash

#SBATCH --job-name=RunFixCluster
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=RunFixCluster.out
#SBATCH --error=RunFixCluster.err
#SBATCH --partition=broadwl
#SBATCH --mem=36G
#SBATCH --mail-type=END


source ~/activate_anaconda.sh
conda activate comp_threeprime_env

gunzip ../data/DiffSplice_liftedJunc/MergeCombined_perind_numers.counts.gz

python fixLeafCluster.py ../data/DiffSplice_liftedJunc/MergeCombined_perind_numers.counts ../data/DiffSplice_liftedJunc/MergeCombined_perind_numers.counts.fixed

gzip ../data/DiffSplice_liftedJunc/MergeCombined_perind_numers.counts.fixed


gunzip ../data/DiffSplice_liftedJunc/MergeCombined_perind.counts.gz

python fixLeafCluster.py ../data/DiffSplice_liftedJunc/MergeCombined_perind.counts ../data/DiffSplice_liftedJunc/MergeCombined_perind.counts.fixed

gzip ../data/DiffSplice_liftedJunc/MergeCombined_perind.counts.fixed
