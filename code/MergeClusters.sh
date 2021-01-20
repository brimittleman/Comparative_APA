#!/bin/bash

#SBATCH --job-name=MergeClusters
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=MergeClusters.out
#SBATCH --error=MergeClusters.err
#SBATCH --partition=broadwl
#SBATCH --mem=16G
#SBATCH --mail-type=END

source ~/activate_anaconda_python2.sh
#conda activate comp_threeprime_env
#source  activate leafcutter


python /project2/yangili1/yangili/leafcutter_scripts/merge_leafcutter_clusters.py ../data/DiffSplice_liftedJunc/CombinedClusters ../data/DiffSplice_liftedJunc/chimpJunc_refined ../data/DiffSplice_liftedJunc/humanJunc_refined
