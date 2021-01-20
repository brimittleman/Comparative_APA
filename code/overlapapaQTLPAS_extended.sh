#!/bin/bash

#SBATCH --job-name=intersectAnnoExt
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=intersectAnnoExt.out
#SBATCH --error=intersectAnnoExt.err
#SBATCH --partition=broadwl
#SBATCH --mem=36G
#SBATCH --mail-type=END

source ~/activate_anaconda.sh
conda activate comp_threeprime_env

bedtools intersect -sorted  -S -wa -a ../data/CompapaQTLpas/CompAPA_PAS_5percHuman.sort.bed -b ../data/liftover_files/APAPAS_GeneLocAnno.5perc.hg19lifted_extended.sort.bed   > ../data/CompapaQTLpas/PAS_5percHuman.sort.Intersect_ext.bed


bedtools intersect -sorted  -v -S -wa -a ../data/CompapaQTLpas/CompAPA_PAS_5percHuman.sort.bed -b ../data/liftover_files/APAPAS_GeneLocAnno.5perc.hg19lifted_extended.sort.bed   > ../data/CompapaQTLpas/PAS_5percHuman.sort.Intersect.NoOverlap_ext.bed
