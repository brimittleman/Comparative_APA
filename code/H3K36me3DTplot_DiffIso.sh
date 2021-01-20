#!/bin/bash

#SBATCH --job-name=H3K36me3DTplot_DiffIso.sh
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=H3K36me3DTplot_DiffIso.out
#SBATCH --error=H3K36me3DTplot_DiffIso.err
#SBATCH --partition=bigmem2
#SBATCH --mem=64G
#SBATCH --mail-type=END


source ~/activate_anaconda.sh
conda activate comp_threeprime_env


computeMatrix scale-regions -S ../data/H3K36me3/MergedIndiv.H3k36me3.sort.bw   -R ../data/PAS_doubleFilter/PAS_diffUsed.bed -b 1000  -a 1000 --skipZeros --transcript_id_designator 4 -out ../data/H3K36me3/H3k36me3PASDiffUsed.gz

plotHeatmap --sortRegions descend -m ../data/H3K36me3/H3k36me3PASDiffUsed.gz --plotTitle "Diff Used PAS H3K36me3" --heatmapHeight 4 --colorMap magma  --averageTypeSummaryPlot "mean" --startLabel "5'"  --endLabel "3'" -out ../data/H3K36me3/H3k36me3PASDiffUsed.png


computeMatrix scale-regions -S ../data/H3K36me3/MergedIndiv.H3k36me3.sort.bw   -R ../data/PAS_doubleFilter/PAS_NOTdiffUsed.bed -b 1000  -a 1000 --skipZeros --transcript_id_designator 4 -out ../data/H3K36me3/H3k36me3PASNOTDiffUsed.gz

plotHeatmap --sortRegions descend -m ../data/H3K36me3/H3k36me3PASNOTDiffUsed.gz --plotTitle "Conserved PAS H3K36me3" --heatmapHeight 4 --colorMap magma  --averageTypeSummaryPlot "mean" --startLabel "5'"  --endLabel "3'" -out ../data/H3K36me3/H3k36me3PASNOTDiffUsed.png
