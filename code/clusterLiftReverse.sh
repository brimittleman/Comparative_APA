#!/bin/bash

#SBATCH --job-name=LiftClustersSecond
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=LiftClustersSecond.out
#SBATCH --error=LiftClustersSecond.err
#SBATCH --partition=broadwl
#SBATCH --mem=36G
#SBATCH --mail-type=END

source ~/activate_anaconda.sh
conda activate comp_threeprime_env


#liftover  oldFile map.chain newFile unMapped


#human to chimp

liftOver ../Chimp/data/RNAseq/DiffSplice/chimpJunc_inHuman.bed ../data/chainFiles/hg38ToPanTro5.over.chain  ../Chimp/data/RNAseq/DiffSplice/chimpJunc_inHuman_B2Chimp.bed  ../Chimp/data/RNAseq/DiffSplice/chimpJunc_inHuman_B2Chimp_unlifted.bed


#chimp to human

liftOver ../Human/data/RNAseq/DiffSplice/humanJunc_inChimp.bed ../data/chainFiles/panTro5ToHg38.over.chain ../Human/data/RNAseq/DiffSplice/humanJunc_inChimp_B2Human.bed ../Human/data/RNAseq/DiffSplice/humanJunc_inChimp_B2Human_unlifted.bed
