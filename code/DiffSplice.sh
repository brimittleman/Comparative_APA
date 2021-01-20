#!/bin/bash

#SBATCH --job-name=DiffSplice
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=DiffSplice.out
#SBATCH --error=DiffSplice.err
#SBATCH --partition=broadwl
#SBATCH --mem=36G
#SBATCH --mail-type=END

module load Anaconda3/5.3.0
source  activate leafcutter
#source ~/activate_anaconda_python2.sh

Rscript ../../leafcutter/scripts/leafcutter_ds.R --num_threads 4 --exon_file=../data/DiffSplice/hg38_ncbiRefseq_exonsfixed.gz ../data/DiffSplice_liftedJunc/MergeCombined_perind_numers.counts.fixed.gz ../data/DiffSplice_liftedJunc/groups_file.txt -o ../data/DiffSplice_liftedJunc/MergedRes
