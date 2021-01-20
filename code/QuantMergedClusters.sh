#!/bin/bash

#SBATCH --job-name=QuantMergeClusters
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=QuantMergeClusters.out
#SBATCH --error=QuantMergeClusters.err
#SBATCH --partition=broadwl
#SBATCH --mem=16G
#SBATCH --mail-type=END

source ~/activate_anaconda_python2.sh
#conda activate comp_threeprime_env
#source  activate leafcutter


python /project2/yangili1/yangili/leafcutter_scripts/leafcutter_merge_regtools.py -m 10 -j ../data/DiffSplice_liftedJunc/BothSpec_juncfiles.txt -o MergeCombined -r ../data/DiffSplice_liftedJunc/ -c ../data/DiffSplice_liftedJunc/CombinedClusters
