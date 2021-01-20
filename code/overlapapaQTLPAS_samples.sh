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

bedtools intersect -sorted  -S -wa -a ../data/OverlapBenchmark/sample1PAS_sort.bed -b ../data/liftover_files/APAPAS_GeneLocAnno.5perc.hg19lifted_extended.sort.bed    > ../data/OverlapBenchmark/sample1PAS_sort.Intersect.bed


bedtools intersect -sorted  -v -S -wa -a ../data/OverlapBenchmark/sample1PAS_sort.bed -b ../data/liftover_files/APAPAS_GeneLocAnno.5perc.hg19lifted_extended.sort.bed    > ../data/OverlapBenchmark/sample1PAS_sort.Intersect.NoOverlap.bed

bedtools intersect -sorted  -S -wa -a ../data/OverlapBenchmark/sample2PAS_sort.bed -b ../data/liftover_files/APAPAS_GeneLocAnno.5perc.hg19lifted_extended.sort.bed    > ../data/OverlapBenchmark/sample2PAS_sort.Intersect.bed


bedtools intersect -sorted  -v -S -wa -a ../data/OverlapBenchmark/sample2PAS_sort.bed -b ../data/liftover_files/APAPAS_GeneLocAnno.5perc.hg19lifted_extended.sort.bed    > ../data/OverlapBenchmark/sample2PAS_sort.Intersect.NoOverlap.bed

bedtools intersect -sorted  -S -wa -a ../data/OverlapBenchmark/sample3PAS_sort.bed -b ../data/liftover_files/APAPAS_GeneLocAnno.5perc.hg19lifted_extended.sort.bed    > ../data/OverlapBenchmark/sample3PAS_sort.Intersect.bed


bedtools intersect -sorted  -v -S -wa -a ../data/OverlapBenchmark/sample3PAS_sort.bed -b ../data/liftover_files/APAPAS_GeneLocAnno.5perc.hg19lifted_extended.sort.bed    > ../data/OverlapBenchmark/sample3PAS_sort.Intersect.NoOverlap.bed
