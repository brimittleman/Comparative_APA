#!/bin/bash

#SBATCH --job-name=LiftClustersFirst
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=LexiftClustersFirst.out
#SBATCH --error=LiftClustersFirst.err
#SBATCH --partition=broadwl
#SBATCH --mem=36G
#SBATCH --mail-type=END

source ~/activate_anaconda.sh
conda activate comp_threeprime_env


#liftover  oldFile map.chain newFile unMapped


#human to chimp

liftOver ../Human/data/RNAseq/DiffSplice/humanJuncNamed.bed ../data/chainFiles/hg38ToPanTro5.over.chain  ../Human/data/RNAseq/DiffSplice/humanJunc_inChimp.bed  ../Human/data/RNAseq/DiffSplice/humanJunc_unlifted.bed


#chimp to human

liftOver ../Chimp/data/RNAseq/DiffSplice/chimpJuncNamed.bed ../data/chainFiles/panTro5ToHg38.over.chain ../Chimp/data/RNAseq/DiffSplice/chimpJunc_inHuman.bed ../Chimp/data/RNAseq/DiffSplice/chimpJunc_unlifted.bed
