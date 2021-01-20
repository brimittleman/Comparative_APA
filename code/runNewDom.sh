#!/bin/bash


#SBATCH --job-name=RunNewDom
#SBATCH --output=RunNewDom.out
#SBATCH --error=RunNewDom.err
#SBATCH --time=1:00:00
#SBATCH --partition=broadwl
#SBATCH --mem=8G
#SBATCH --mail-type=END


source ~/activate_anaconda.sh
conda activate comp_threeprime_env

for i in .1 .2 .3 .4 .5 .6 .7 .8 .9
do
echo $i
python FindDomXCutoff.py ../data/PAS_doubleFilter/PAS_5perc_either_HumanCoord_BothUsage_meta_doubleFilter.txt ../data/DomDefGreaterX/Human_${i}_dominantPAS.txt $i Human
done

for i in .1 .2 .3 .4 .5 .6 .7 .8 .9
do
echo $i
python FindDomXCutoff.py ../data/PAS_doubleFilter/PAS_5perc_either_HumanCoord_BothUsage_meta_doubleFilter.txt ../data/DomDefGreaterX/Chimp_${i}_dominantPAS.txt $i Chimp
done
