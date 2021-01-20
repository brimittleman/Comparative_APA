#!/bin/bash


#SBATCH --job-name=liftoverPAShg19to38
#SBATCH --time=24:00:00
#SBATCH --partition=broadwl
#SBATCH --mem=16G
#SBATCH --tasks-per-node=4
#SBATCH --mail-type=END
#SBATCH --output=liftoverPAShg19to38.out
#SBATCH --error=liftoverPAShg19to38.err



source ~/activate_anaconda.sh
conda activate comp_threeprime_env

liftOver ../data/liftover_files/APAPAS_withCHR_GeneLocAnno.5perc.sort.bed  ../data/liftover_files/hg19ToHg38.over.chain.gz ../data/liftover_files/APAPAS_GeneLocAnno.5perc.hg19lifted.bed ../data/liftover_files/APAPAS_GeneLocAnno.5perc.hg19unlifted.bed
