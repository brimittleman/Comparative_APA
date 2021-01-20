#!/bin/bash

#SBATCH --job-name=runCountNucleotidesPantro6
#SBATCH --time=24:00:00
#SBATCH --partition=broadwl
#SBATCH --mem=4G
#SBATCH --tasks-per-node=4
#SBATCH --mail-type=END
#SBATCH --output=runCountNucleotidesPantro6.out
#SBATCH --error=runCountNucleotidesPantro6.err



source ~/activate_anaconda.sh
conda activate comp_threeprime_env

touch  ../data/EvalPantro5/CountsforNucleotides_pantro6.txt


for i in  $(ls /project2/gilad/briana/CompAPA_pantro6/Human/data/nuc_10up/*)
do
python CountNucleotides.py $i "N" >> ../data/EvalPantro5/CountsforNucleotides_pantro6.txt
done

for i in $(ls /project2/gilad/briana/CompAPA_pantro6/Human/data/nuc_10up/*)
do
python CountNucleotides.py $i "T" >> ../data/EvalPantro5/CountsforNucleotides_pantro6.txt
done

for i in $(ls /project2/gilad/briana/CompAPA_pantro6/Chimp/data/nuc_10up/*)
do
python CountNucleotides.py $i "N" >> ../data/EvalPantro5/CountsforNucleotides_pantro6.txt
done

for i in $(ls /project2/gilad/briana/CompAPA_pantro6/Chimp/data/nuc_10up/*)
do
python CountNucleotides.py $i "T" >> ../data/EvalPantro5/CountsforNucleotides_pantro6.txt
done
