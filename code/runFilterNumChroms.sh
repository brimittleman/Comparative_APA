#!/bin/bash

#SBATCH --job-name=filterJuncChroms
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=filterJuncChroms.out
#SBATCH --error=filterJuncChroms.err
#SBATCH --partition=broadwl
#SBATCH --mem=16G
#SBATCH --mail-type=END

source ~/activate_anaconda.sh
conda activate comp_threeprime_env

touch ../data/DiffSplice_liftedJunc/BothSpec_juncfiles.txt

for juncfile in `ls /project2/gilad/briana/Comparative_APA/Human/data/RNAseq/sort/*.backHuman.junc.SamePlace`
do
    echo $juncfile
    python filterNumChroms.py $juncfile
    echo $juncfile.fixed.junc >> ../data/DiffSplice_liftedJunc/BothSpec_juncfiles.txt

done



for juncfile in `ls /project2/gilad/briana/Comparative_APA/Chimp/data/RNAseq/sort/*SamePlace.FinalInHuman.junc`
do
    echo $juncfile
    python filterNumChroms.py $juncfile
    echo $juncfile.fixed.junc >> ../data/DiffSplice_liftedJunc/BothSpec_juncfiles.txt
done
