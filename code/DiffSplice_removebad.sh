#!/bin/bash

#SBATCH --job-name=DiffSplice_removebad
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=DiffSplice_removebad.out
#SBATCH --error=DiffSplice_removebad.err
#SBATCH --partition=broadwl
#SBATCH --mem=36G
#SBATCH --mail-type=END

source deactivate
source deactivate
module load Anaconda3/5.3.0
source  activate leafcutter




 Rscript ../../leafcutter/scripts/leafcutter_ds.R --num_threads 4 --exon_file=../data/DiffSplice/hg38_ncbiRefseq_exonsfixed.gz ../data/DiffSplice_removeBad/BothSpec_perind.counts.gz ../data/DiffSplice_removeBad/groups_file.txt -o ../data/DiffSplice_removeBad/
