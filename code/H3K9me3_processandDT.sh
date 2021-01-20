#!/bin/bash

#SBATCH --job-name=H3K9me3_processandDT
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=mergeandsort_H3K9me3
#SBATCH --error=mergeandsort_H3K9me3
#SBATCH --partition=bigmem2
#SBATCH --mem=100G
#SBATCH --mail-type=END

source ~/activate_anaconda.sh
conda activate comp_threeprime_env

samtools merge ../data/H3K36me3/MergedIndiv.H3K9me3.bam  ../data/H3K36me3/*H3K9me3.bam

samtools sort ../data/H3K36me3/MergedIndiv.H3K9me3.bam  -o ../data/H3K36me3/MergedIndiv.H3K9me3.sort.bam 

samtools index ../data/H3K36me3/MergedIndiv.H3K9me3.sort.bam 

bamCoverage -b ../data/H3K36me3/MergedIndiv.H3K9me3.sort.bam  -o ../data/H3K36me3/MergedIndiv.H3K9me3.sort.bw  


computeMatrix scale-regions -S ../data/H3K36me3/MergedIndiv.H3K9me3.sort.bw   -R ../data/PAS_doubleFilter/PAS_doublefilter_either_HumanCoordHummanUsage.sort.bed -b 500  -a 500 --skipZeros --transcript_id_designator 4 -out ../data/H3K36me3/H3K9me3PAS.gz

plotHeatmap --sortRegions descend -m ../data/H3K36me3/H3K9me3PAS.gz --plotTitle "PAS H3K9me3" --heatmapHeight 4 --colorMap magma  --averageTypeSummaryPlot "mean" --startLabel "5'"  --endLabel "3'" -out ../data/H3K36me3/H3K9me3PAS.png


computeMatrix scale-regions -S ../data/H3K36me3/MergedIndiv.H3K9me3.sort.bw   -R ../data/PAS_doubleFilter/PAS_doublefilter_either_HumanCoordHummanUsage.sort.bed -b 2000  -a 2000 --skipZeros --transcript_id_designator 4 -out ../data/H3K36me3/H3K9me3PASwide.gz

plotHeatmap --sortRegions descend -m ../data/H3K36me3/H3K9me3PASwide.gz --plotTitle "PAS H3K9me3" --heatmapHeight 4 --colorMap magma  --averageTypeSummaryPlot "mean" --startLabel "5'"  --endLabel "3'" -out ../data/H3K36me3/H3K9me3PASwide.png

computeMatrix scale-regions -S ../data/H3K36me3/MergedIndiv.H3K9me3.sort.bw   -R ../data/PAS_doubleFilter/DistalPAS_orthoUTR.bed -b 10000  -a 10000 --skipZeros --transcript_id_designator 4 -out ../data/H3K36me3/H3K9me3distalPAS_longer.gz

plotHeatmap --sortRegions descend -m ../data/H3K36me3/H3K9me3distalPAS_longer.gz --plotTitle "Distal UTR PAS H3K9me3" --heatmapHeight 4 --colorMap magma  --averageTypeSummaryPlot "mean" --startLabel "5'"  --endLabel "3'" -out ../data/H3K36me3/H3K9me3distalPAS_longer.png


computeMatrix scale-regions -S ../data/H3K36me3/MergedIndiv.H3K9me3.sort.bw   -R ../data/PAS_doubleFilter/PAS_doublefilter_either_HumanCoordHummanUsage.sort.bed -b 2000  -a 2000 --skipZeros --transcript_id_designator 4 -out ../data/H3K36me3/H3K9me3PASwide.gz

plotHeatmap --sortRegions descend -m ../data/H3K36me3/H3K9me3PASwide.gz --plotTitle "PAS H3K9me3" --heatmapHeight 4 --colorMap magma  --averageTypeSummaryPlot "mean" --startLabel "5'"  --endLabel "3'" -out ../data/H3K36me3/H3K9me3PASwide.png
