#!/bin/bash

#SBATCH --job-name=JunctionLiftFinalChimp
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=JunctionLiftFinalChimp.out
#SBATCH --error=JunctionLiftFinalChimp.err
#SBATCH --partition=broadwl
#SBATCH --mem=36G
#SBATCH --mail-type=END

source ~/activate_anaconda.sh
conda activate comp_threeprime_env



for juncfile in `ls /project2/gilad/briana/Comparative_APA/Chimp/data/RNAseq/sort/*.junc.2Human.backChimp.SamePlace`
do
    echo $juncfile
    liftOver $juncfile ../data/chainFiles/panTro6ToHg38.over.chain  $juncfile.FinalInHuman.junc $juncfile.unlift -bedPlus=12
done
