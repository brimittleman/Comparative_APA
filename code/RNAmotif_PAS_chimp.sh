#!/bin/bash

#SBATCH --job-name=RNAmotif_PAS_chimp
#SBATCH --account=pi-gilad
#SBATCH --time=36:00:00
#SBATCH --output=RNAmotif_PAS_chimp.out
#SBATCH --error=RNAmotif_PAS_chimp.err
#SBATCH --partition=bigmem2
#SBATCH --mem=100G
#SBATCH --mail-type=END


source ~/activate_anaconda.sh
conda activate comp_threeprime_env


findMotifsGenome.pl ../data/PAS_doubleFilter/PAS_doublefilter_either_ChimpCoordChimpUsage.sort.bed /project2/gilad/briana/genome_anotation_data/Chimp_genome/panTro6.fa  ../data/PAS_doubleFilter/FindMotif_chimp/  -rna -h -len 6
