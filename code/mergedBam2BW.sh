#!/bin/bash

#SBATCH --job-name=mergedbam2bw
#SBATCH --account=pi-yangili1
#SBATCH --time=36:00:00
#SBATCH --output=mergedbam2bw.out
#SBATCH --error=mergedbam2bw.err
#SBATCH --cpus-per-task=2
#SBATCH --partition=bigmem2
#SBATCH --mem=216G
#SBATCH --mail-type=END

source ~/activate_anaconda.sh
conda activate comp_threeprime_env




bamCoverage -b ../Human/data/mergedbyFracBam/human_Nuclear.SamplesMerged.sort.bam  -o  ../Human/data/mergedbw_byFrac/human_Nuclear.SamplesMerged.sort.bw

bamCoverage -b ../Human/data/mergedbyFracBam/human_Total.SamplesMerged.sort.bam  -o  ../Human/data/mergedbw_byFrac/human_Total.SamplesMerged.sort.bw

bamCoverage -b ../Chimp/data/mergedbyFracBam/chimp_Nuclear.SamplesMerged.sort.bam  -o  ../Chimp/data/mergedbw_byFrac/chimp_Nuclear.SamplesMerged.sort.bw

bamCoverage -b ../Chimp/data/mergedbyFracBam/chimp_Total.SamplesMerged.sort.bam  -o  ../Chimp/data/mergedbw_byFrac/chimp_Total.SamplesMerged.sort.bw
