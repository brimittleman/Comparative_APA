#!/bin/bash

#SBATCH --job-name=RunPosMCMediationDF
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=RunPosMCMediationDF.out
#SBATCH --error=RunPosMCMediationDF.err
#SBATCH --partition=broadwl
#SBATCH --mem=36G
#SBATCH --mail-type=END


source ~/activate_anaconda.sh
conda activate comp_threeprime_env


Rscript postiveMediation_montecarlo_DF.R
