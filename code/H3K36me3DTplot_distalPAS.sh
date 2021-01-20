#!/bin/bash

#SBATCH --job-name=H3K36me3DTplot_distalPAS.sh
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=H3K36me3DTplot_distalPAS.out
#SBATCH --error=H3K36me3DTplot_distalPAS.err
#SBATCH --partition=bigmem2
#SBATCH --mem=100G
#SBATCH --mail-type=END


source ~/activate_anaconda.sh
conda activate comp_threeprime_env


computeMatrix scale-regions -S ../data/H3K36me3/MergedIndiv.H3k36me3.sort.bw   -R ../data/PAS_doubleFilter/DistalPAS_orthoUTR.bed -b 2000  -a 2000 --skipZeros --transcript_id_designator 4 -out ../data/H3K36me3/H3k36me3distalPAS.gz

plotHeatmap --sortRegions descend -m ../data/H3K36me3/H3k36me3distalPAS.gz --plotTitle "Distal UTR PAS H3L36me3" --heatmapHeight 4 --colorMap magma  --averageTypeSummaryPlot "mean" --startLabel "5'"  --endLabel "3'" -out ../data/H3K36me3/H3k36me3distalPAS.png


computeMatrix scale-regions -S ../data/H3K36me3/MergedIndiv.H3k36me3.sort.bw   -R ../data/PAS_doubleFilter/DistalPAS_orthoUTR.bed -b 10000  -a 10000 --skipZeros --transcript_id_designator 4 -out ../data/H3K36me3/H3k36me3distalPAS_longer.gz

plotHeatmap --sortRegions descend -m ../data/H3K36me3/H3k36me3distalPAS_longer.gz --plotTitle "Distal UTR PAS H3L36me3" --heatmapHeight 4 --colorMap magma  --averageTypeSummaryPlot "mean" --startLabel "5'"  --endLabel "3'" -out ../data/H3K36me3/H3k36me3distalPAS_longer.png


