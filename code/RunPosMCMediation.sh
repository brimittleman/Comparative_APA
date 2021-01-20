#!/bin/bash

#SBATCH --job-name=RunPosMCMediation
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=RunPosMCMediationr.out
#SBATCH --error=RunPosMCMediation.err
#SBATCH --partition=broadwl
#SBATCH --mem=36G
#SBATCH --mail-type=END


source ~/activate_anaconda.sh
conda activate comp_threeprime_env


Rscript postiveMediation_montecarlo.R
