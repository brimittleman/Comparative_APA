#!/bin/bash

#SBATCH --job-name=runCountNucleotides
#SBATCH --time=24:00:00
#SBATCH --partition=broadwl
#SBATCH --mem=4G
#SBATCH --tasks-per-node=4
#SBATCH --mail-type=END
#SBATCH --output=runCountNucleotides.out
#SBATCH --error=runCountNucleotides.err



source ~/activate_anaconda.sh
conda activate comp_threeprime_env

touch  ../data/EvalPantro5/CountsforNucleotides.txt


for i in  $(ls /project2/gilad/briana/Comparative_APA/Human/data/nuc_10up/*)
do
python CountNucleotides.py $i "N" >> ../data/EvalPantro5/CountsforNucleotides.txt
done

for i in $(ls /project2/gilad/briana/Comparative_APA/Human/data/nuc_10up/*)
do
python CountNucleotides.py $i "T" >> ../data/EvalPantro5/CountsforNucleotides.txt
done

for i in $(ls /project2/gilad/briana/Comparative_APA/Chimp/data/nuc_10up/*)
do
python CountNucleotides.py $i "N" >> ../data/EvalPantro5/CountsforNucleotides.txt
done

for i in $(ls /project2/gilad/briana/Comparative_APA/Chimp/data/nuc_10up/*)
do
python CountNucleotides.py $i "T" >> ../data/EvalPantro5/CountsforNucleotides.txt
done
