#!/bin/bash

#SBATCH --job-name=RunNegMCMediation
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=RunNegMCMediationr.out
#SBATCH --error=RunNegMCMediation.err
#SBATCH --partition=broadwl
#SBATCH --mem=36G
#SBATCH --mail-type=END


source ~/activate_anaconda.sh
conda activate comp_threeprime_env


Rscript negativeMediation_montecarlo.R
