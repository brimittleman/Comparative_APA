#!/bin/bash

#SBATCH --job-name=quatJunc
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=quatJunc.out
#SBATCH --error=quatJunc.err
#SBATCH --partition=broadwl
#SBATCH --mem=16G
#SBATCH --mail-type=END

source  ~/activate_anaconda_python2.sh
#source  activate leafcutter



python /project2/gilad/briana/leafcutter/clustering/leafcutter_cluster.py -j /project2/gilad/briana/Comparative_APA/Human/data/RNAseq/DiffSplice_removeBad/human_juncfiles.txt -m 50 -o humanJunc -l 500000

mv *humanJunc* /project2/gilad/briana/Comparative_APA/Human/data/RNAseq/DiffSplice_removeBad/


python /project2/gilad/briana/leafcutter/clustering/leafcutter_cluster.py -j /project2/gilad/briana/Comparative_APA/Chimp/data/RNAseq/DiffSplice_removeBad/chimp_juncfiles.txt -m 50 -o chimpJunc -l 500000

mv *chimpJunc* /project2/gilad/briana/Comparative_APA/Chimp/data/RNAseq/DiffSplice_removeBad/
