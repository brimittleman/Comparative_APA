#!/bin/bash

#SBATCH --job-name=LiftClustersFirst_remove
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=LiftClustersFirst_remove.out
#SBATCH --error=LiftClustersFirst_remove.err
#SBATCH --partition=broadwl
#SBATCH --mem=36G
#SBATCH --mail-type=END

source ~/activate_anaconda.sh
conda activate comp_threeprime_env


#liftover  oldFile map.chain newFile unMapped


#human to chimp

liftOver ../Human/data/RNAseq/DiffSplice_removeBad/humanJuncNamed.bed ../data/chainFiles/hg38ToPanTro5.over.chain  ../Human/data/RNAseq/DiffSplice_removeBad/humanJunc_inChimp.bed  ../Human/data/RNAseq/DiffSplice_removeBad/humanJunc_unlifted.bed


#chimp to human

liftOver ../Chimp/data/RNAseq/DiffSplice_removeBad/chimpJuncNamed.bed ../data/chainFiles/panTro5ToHg38.over.chain ../Chimp/data/RNAseq/DiffSplice_removeBad/chimpJunc_inHuman.bed ../Chimp/data/RNAseq/DiffSplice_removeBad/chimpJunc_unlifted.bed
