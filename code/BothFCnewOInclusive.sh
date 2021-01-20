#!/bin/bash

#SBATCH --job-name=BothFCnewOInclusive
#SBATCH --time=24:00:00
#SBATCH --partition=broadwl
#SBATCH --mem=32G
#SBATCH --mail-type=END
#SBATCH --output=BothFCnewOInclusive.sh.out
#SBATCH --error=BothFCnewOInclusive.sh.err



source ~/activate_anaconda.sh
conda activate comp_threeprime_env
#actual strand like in original

#featureCounts -a {input.annotation} -F SAF -o {output} {input.inputdir}*.sort.bam -s 1

featureCounts -a ../Human/data/inclusivePeaks/human_APApeaks.ALLChrom.SAF  -F SAF -o ../data/testQuant/Human_AllPeaksInclusive_newFC  ../Human/data/sort_clean/*_N-clean.sort.bam -s 1 -O

featureCounts -a ../Chimp/data/inclusivePeaks/chimp_APApeaks.ALLChrom.SAF  -F SAF -o ../data/testQuant/Chimp_AllPeaksInclusive_newFC  ../Chimp/data/sort_clean/*_N-clean.sort.bam -s 1 -O
