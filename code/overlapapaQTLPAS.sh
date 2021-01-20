#!/bin/bash

#SBATCH --job-name=intersectAnno
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=intersectAnno.out
#SBATCH --error=intersectAnno.err
#SBATCH --partition=broadwl
#SBATCH --mem=36G
#SBATCH --mail-type=END

source ~/activate_anaconda.sh
conda activate comp_threeprime_env

bedtools intersect -sorted  -S -wa -a ../data/CompapaQTLpas/CompAPA_PAS_5percHuman.sort.bed -b ../data/liftover_files/APAPAS_GeneLocAnno.5perc.hg19lifted.sorted.bed   > ../data/CompapaQTLpas/PAS_5percHuman.sort.Intersect.bed


bedtools intersect -sorted  -v -S -wa -a ../data/CompapaQTLpas/CompAPA_PAS_5percHuman.sort.bed -b ../data/liftover_files/APAPAS_GeneLocAnno.5perc.hg19lifted.sorted.bed   > ../data/CompapaQTLpas/PAS_5percHuman.sort.Intersect.NoOverlap.bed
