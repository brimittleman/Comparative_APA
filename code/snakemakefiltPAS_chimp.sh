#!/bin/bash

# sbatch submission script to run main snakemake process. It then submits
# individual jobs from the compute node.


#SBATCH --job-name=snakemakePASFiltchimp
#SBATCH --time=24:00:00
#SBATCH --partition=broadwl
#SBATCH --mem=4G
#SBATCH --tasks-per-node=4
#SBATCH --mail-type=END
#SBATCH --output=snakemakePASFiltChimp.out
#SBATCH --error=snakemakePASFiltChimp.err



source ~/activate_anaconda.sh
conda activate comp_threeprime_env


bash submit-snakemakefiltPAS-chimp.sh $*
