#!/bin/bash

#SBATCH --job-name=quatJunc
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=quatJunc.out
#SBATCH --error=quatJunc.err
#SBATCH --partition=broadwl
#SBATCH --mem=16G
#SBATCH --mail-type=END

source ~/activate_anaconda_python2.sh
#conda activate comp_threeprime_env
#source  activate leafcutter



python /project2/gilad/briana/leafcutter/clustering/leafcutter_cluster_regtools.py -j ../data/DiffSplice_liftedJunc/Human_juncfiles.txt -m 10 -o humanJunc -l 500000

mv *humanJunc* /project2/gilad/briana/Comparative_APA/data/DiffSplice_liftedJunc/

python /project2/gilad/briana/leafcutter/clustering/leafcutter_cluster_regtools.py -j ../data/DiffSplice_liftedJunc/Chimp_juncfiles.txt -m 10 -o chimpJunc -l 500000


mv *chimpJunc* /project2/gilad/briana/Comparative_APA/data/DiffSplice_liftedJunc/
