#!/bin/bash

#SBATCH --job-name=BothFCMM
#SBATCH --time=24:00:00
#SBATCH --partition=broadwl
#SBATCH --mem=32G
#SBATCH --mail-type=END
#SBATCH --output=BothFCMM.out
#SBATCH --error=BothFCMM.err



source ~/activate_anaconda.sh
conda activate comp_threeprime_env


#not double filtered


featureCounts -O -M -a ../data/cleanPeaks_anno/AllPAS_postLift.sort_LocAnnoPARSED_chimpLoc.SAF -F SAF -o ../data/testQuant/ALLPAS_postLift_LocParsed_Chimp_Mflag ../Chimp/data/sort_clean/*_N-clean.sort.bam  -s 1

featureCounts -O --primary -a ../data/cleanPeaks_anno/AllPAS_postLift.sort_LocAnnoPARSED_chimpLoc.SAF -F SAF -o ../data/testQuant/ALLPAS_postLift_LocParsed_Chimp_Primary ../Chimp/data/sort_clean/*_N-clean.sort.bam  -s 1


featureCounts -O -M -a ../data/cleanPeaks_anno/AllPAS_postLift.sort_LocAnnoPARSED_chimpLoc.SAF -F SAF -o ../data/testQuant/ALLPAS_postLift_LocParsed_Chimp_BothFrac_Mflag ../Chimp/data/sort_clean/*.bam  -s 1

featureCounts -O -M -a ../data/cleanPeaks_anno/AllPAS_postLift.sort_LocAnnoPARSED.SAF -F SAF -o ../data/testQuant/ALLPAS_postLift_LocParsed_Human_Mflag ../Human/data/sort_clean/*_N-clean.sort.bam -s 1

featureCounts -O -M -a ../data/cleanPeaks_anno/AllPAS_postLift.sort_LocAnnoPARSED.SAF -F SAF -o ../data/testQuant/ALLPAS_postLift_LocParsed_Human_BothFrac_Mflag ../Human/data/sort_clean/*.bam -s 1

featureCounts -O --primary -a ../data/cleanPeaks_anno/AllPAS_postLift.sort_LocAnnoPARSED.SAF -F SAF -o ../data/testQuant/ALLPAS_postLift_LocParsed_Human_Primary ../Human/data/sort_clean/*_N-clean.sort.bam -s 1
