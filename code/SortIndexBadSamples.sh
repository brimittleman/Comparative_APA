#!/bin/bash

#SBATCH --job-name=SortIndexBadSamples
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=SortIndexBadSamples.out
#SBATCH --error=SortIndexBadSamples.err
#SBATCH --partition=broadwl
#SBATCH --mem=36G
#SBATCH --mail-type=END

source ~/activate_anaconda.sh
conda activate comp_threeprime_env

samtools sort ../data/TwoBadSampleAnalysis/NA4973inHuman > ../data/TwoBadSampleAnalysis/NA4973inHuman-sort.bam

samtools index ../data/TwoBadSampleAnalysis/NA4973inHuman-sort.bam

samtools sort ../data/TwoBadSampleAnalysis/NA18498inChimp > ../data/TwoBadSampleAnalysis/NA18498inChimp-sort.bam

samtools index ../data/TwoBadSampleAnalysis/NA18498inChimp-sort.bam
