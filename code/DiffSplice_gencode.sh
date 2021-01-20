#!/bin/bash

#SBATCH --job-name=GencodeDiffSplice
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=GencodeDiffSplice.out
#SBATCH --error=GencodeDiffSplice.err
#SBATCH --partition=broadwl
#SBATCH --mem=36G
#SBATCH --mail-type=END

module load Anaconda3/5.3.0
source  activate leafcutter




 Rscript ../../leafcutter/scripts/leafcutter_ds.R --num_threads 4 --exon_file=../data/DiffSplice/gencode19_exons.txt.gz ../data/DiffSplice/BothSpec_perind.counts.gz ../data/DiffSplice/groups_file.txt -o ../data/DiffSplice/Gencode_
