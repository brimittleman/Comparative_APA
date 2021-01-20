#!/bin/bash


#SBATCH --job-name=wrap_verifybam
#SBATCH --output=wrap_verifybam.out
#SBATCH --error=wrap_verifybam.err
#SBATCH --time=36:00:00
#SBATCH --partition=broadwl
#SBATCH --ntasks-per-node=13
#SBATCH --mem=16G
#SBATCH --mail-type=END


source ~/activate_anaconda.sh
conda activate comp_threeprime_env


for i in $(ls /project2/gilad/briana/Comparative_APA/Human/data/sort_clean/*.bam)
do
sbatch run_verifyBam.sh ${i}
done
