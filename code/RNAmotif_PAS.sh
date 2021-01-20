#!/bin/bash

#SBATCH --job-name=RNAmotif_PAS
#SBATCH --account=pi-gilad
#SBATCH --time=36:00:00
#SBATCH --output=RNAmotif_PAS.out
#SBATCH --error=RNAmotif_PAS.err
#SBATCH --partition=bigmem2
#SBATCH --mem=100G
#SBATCH --mail-type=END


source ~/activate_anaconda.sh
conda activate comp_threeprime_env


findMotifsGenome.pl ../data/PAS_doubleFilter/PAS_doublefilter_either_HumanCoordHummanUsage.bed /project2/gilad/kenneth/References/human/genome/hg38.fa  ../data/PAS_doubleFilter/FindMotif/  -rna -h -len 6
