#!/bin/bash

#SBATCH --job-name=DiffDom_RNAmotif_4
#SBATCH --account=pi-gilad
#SBATCH --time=36:00:00
#SBATCH --output=DiffDom_RNAmotif_4.out
#SBATCH --error=DiffDom_RNAmotif_4.err
#SBATCH --partition=bigmem2
#SBATCH --mem=100G
#SBATCH --mail-type=END


source ~/activate_anaconda.sh
conda activate comp_threeprime_env


findMotifsGenome.pl ../data/DistTwoDom/SeqBetweenDom_4.bed /project2/gilad/kenneth/References/human/genome/hg38.fa  ../data/DistTwoDom/FindMotif  -rna -h -len 8


findMotifs.pl ../data/DistTwoDom/SeqBetweenDom_4_genes.txt human-mRNA ../data/DistTwoDom/mRNAMotif/ -rna -len 8
