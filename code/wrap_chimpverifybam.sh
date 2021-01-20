#!/bin/bash


#SBATCH --job-name=wrap_Chimpverifybam
#SBATCH --output=wrap_Chimpverifybam.out
#SBATCH --error=wrap_Chimpverifybam.err
#SBATCH --time=36:00:00
#SBATCH --partition=broadwl
#SBATCH --ntasks-per-node=13
#SBATCH --mem=16G
#SBATCH --mail-type=END


source ~/activate_anaconda.sh
conda activate comp_threeprime_env


for i in $(ls /project2/gilad/briana/Comparative_APA/Chimp/data/sort_clean/*.bam)
do
sbatch run_chimpverifybam.sh ${i}
done
