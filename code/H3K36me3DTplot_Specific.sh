#!/bin/bash

#SBATCH --job-name=H3K36me3DTplot_Specific.sh
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=H3K36me3DTplot_Specific.out
#SBATCH --error=H3K36me3DTplot_Specific.err
#SBATCH --partition=bigmem2
#SBATCH --mem=64G
#SBATCH --mail-type=END


source ~/activate_anaconda.sh
conda activate comp_threeprime_env


computeMatrix scale-regions -S ../data/H3K36me3/MergedIndiv.H3k36me3.sort.bw -R ../data/PAS_doubleFilter/PAS_HumanSpecifc.bed -b 1000  -a 1000 --skipZeros --transcript_id_designator 4 -out ../data/H3K36me3/H3k36me3HumanSpec.gz

plotHeatmap --sortRegions descend -m ../data/H3K36me3/H3k36me3HumanSpec.gz --plotTitle "Human Specific H3K36me3" --heatmapHeight 4 --colorMap magma  --averageTypeSummaryPlot "mean" --startLabel "5'"  --endLabel "3'" -out ../data/H3K36me3/H3k36me3HumanSpec.png


computeMatrix scale-regions -S ../data/H3K36me3/MergedIndiv.H3k36me3.sort.bw  -R ../data/PAS_doubleFilter/PAS_ChimpSpecifc.bed -b 1000  -a 1000 --skipZeros --transcript_id_designator 4 -out ../data/H3K36me3/H3k36me3ChimpSpec.gz

plotHeatmap --sortRegions descend -m ../data/H3K36me3/H3k36me3ChimpSpec.gz --plotTitle "Chimp Specific H3K36me3" --heatmapHeight 4 --colorMap magma  --averageTypeSummaryPlot "mean" --startLabel "5'"  --endLabel "3'" -out ../data/H3K36me3/H3k36me3ChimpSpec.png
