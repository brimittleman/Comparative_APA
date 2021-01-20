#!/bin/bash

#SBATCH --job-name=DiffDom_RNAmotif_4_splitDE
#SBATCH --account=pi-gilad
#SBATCH --time=36:00:00
#SBATCH --output=DiffDom_RNAmotif_4_splitDE.out
#SBATCH --error=DiffDom_RNAmotif_4_splitDE.err
#SBATCH --partition=bigmem2
#SBATCH --mem=100G
#SBATCH --mail-type=END


source ~/activate_anaconda.sh
conda activate comp_threeprime_env


findMotifsGenome.pl ../data/DistTwoDom/SeqBetweenDom_4_withDE.bed /project2/gilad/kenneth/References/human/genome/hg38.fa  ../data/DistTwoDom/FindMotif/DE_8  -rna -h -len 8 -bg ../data/DistTwoDom/SeqBetweenDom_4_withNODE.bed


findMotifsGenome.pl ../data/DistTwoDom/SeqBetweenDom_4_withDE.bed /project2/gilad/kenneth/References/human/genome/hg38.fa  ../data/DistTwoDom/FindMotif/DE_10  -rna -h -len 10 -bg ../data/DistTwoDom/SeqBetweenDom_4_withNODE.bed


findMotifsGenome.pl ../data/DistTwoDom/SeqBetweenDom_4_withDE.bed /project2/gilad/kenneth/References/human/genome/hg38.fa  ../data/DistTwoDom/FindMotif/DE_12  -rna -h -len 12 -bg ../data/DistTwoDom/SeqBetweenDom_4_withNODE.bed
