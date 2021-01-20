#!/bin/bash

#SBATCH --job-name=RunNegMCMediationDF
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=RunNegMCMediationDF.out
#SBATCH --error=RunNegMCMediationDF.err
#SBATCH --partition=broadwl
#SBATCH --mem=36G
#SBATCH --mail-type=END


source ~/activate_anaconda.sh
conda activate comp_threeprime_env


Rscript negativeMediation_montecarloDF.R
