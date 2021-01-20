#!/bin/bash

#SBATCH --job-name=H3K36me3DTplot.sh
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=H3K36me3DTplot.out
#SBATCH --error=H3K36me3DTplot.err
#SBATCH --partition=bigmem2
#SBATCH --mem=100G
#SBATCH --mail-type=END


source ~/activate_anaconda.sh
conda activate comp_threeprime_env

bamCoverage -b ../data/H3K36me3/MergedIndiv.H3k36me3.sort.bam  -o ../data/H3K36me3/MergedIndiv.H3k36me3.sort.bw  

computeMatrix scale-regions -S ../data/H3K36me3/MergedIndiv.H3k36me3.sort.bw   -R ../data/PAS_doubleFilter/PAS_doublefilter_either_HumanCoordHummanUsage.sort.bed -b 500  -a 500 --skipZeros --transcript_id_designator 4 -out ../data/H3K36me3/H3k36me3PAS.gz

plotHeatmap --sortRegions descend -m ../data/H3K36me3/H3k36me3PAS.gz --plotTitle "PAS H3L36me3" --heatmapHeight 4 --colorMap magma  --averageTypeSummaryPlot "mean" --startLabel "5'"  --endLabel "3'" -out ../data/H3K36me3/H3k36me3PAS.png
