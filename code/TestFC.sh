#!/bin/bash

#SBATCH --job-name=TestFC
#SBATCH --time=24:00:00
#SBATCH --partition=broadwl
#SBATCH --mem=32G
#SBATCH --mail-type=END
#SBATCH --output=TestFC.out
#SBATCH --error=TestFC.err



source ~/activate_anaconda.sh
conda activate comp_threeprime_env


#opposite strand
#../data/PAS_SAF/Chimp_PASdoublefilter.SAF
#../data/PAS_SAF/Human_PASdoublefilter.SAF
#../data/Test_FC_methods
#../data/TestMM2/*MM2.sort.bam


#standard
featureCounts -O -a ../data/PAS_SAF/Chimp_PASdoublefilter.SAF -F SAF -o ../data/Test_FC_methods/Chimp_standard ../data/TestMM2/chimp*N.MM2.sort.bam  -s 2
featureCounts -O -a ../data/PAS_SAF/Human_PASdoublefilter.SAF -F SAF -o ../data/Test_FC_methods/Human_standard ../data/TestMM2/human*N.MM2.sort.bam  -s 2

#multi
featureCounts -O -M -a ../data/PAS_SAF/Chimp_PASdoublefilter.SAF -F SAF -o ../data/Test_FC_methods/Chimp_MM ../data/TestMM2/chimp*N.MM2.sort.bam  -s 2
featureCounts -O -M  -a ../data/PAS_SAF/Human_PASdoublefilter.SAF -F SAF -o ../data/Test_FC_methods/Human_MM ../data/TestMM2/human*N.MM2.sort.bam  -s 2

#primary
featureCounts -O -M --primary -a ../data/PAS_SAF/Chimp_PASdoublefilter.SAF -F SAF -o ../data/Test_FC_methods/Chimp_Primary ../data/TestMM2/chimp*N.MM2.sort.bam  -s 2
featureCounts -O -M --primary   -a ../data/PAS_SAF/Human_PASdoublefilter.SAF -F SAF -o ../data/Test_FC_methods/Human_Primary ../data/TestMM2/human*N.MM2.sort.bam  -s 2

#strigent
featureCounts -O -M -Q 5 -a ../data/PAS_SAF/Chimp_PASdoublefilter.SAF -F SAF -o ../data/Test_FC_methods/Chimp_5 ../data/TestMM2/chimp*N.MM2.sort.bam  -s 2
featureCounts -O -M -Q 5 -a ../data/PAS_SAF/Human_PASdoublefilter.SAF -F SAF -o ../data/Test_FC_methods/Human_5 ../data/TestMM2/human*N.MM2.sort.bam  -s 2


#med
featureCounts -O -M -Q 10 -a ../data/PAS_SAF/Chimp_PASdoublefilter.SAF -F SAF -o ../data/Test_FC_methods/Chimp_10 ../data/TestMM2/chimp*N.MM2.sort.bam  -s 2
featureCounts -O -M -Q 10 -a ../data/PAS_SAF/Human_PASdoublefilter.SAF -F SAF -o ../data/Test_FC_methods/Human_10 ../data/TestMM2/human*N.MM2.sort.bam  -s 2


#lax

featureCounts -O -M -Q 20 -a ../data/PAS_SAF/Chimp_PASdoublefilter.SAF -F SAF -o ../data/Test_FC_methods/Chimp_20 ../data/TestMM2/chimp*N.MM2.sort.bam  -s 2
featureCounts -O -M -Q 20 -a ../data/PAS_SAF/Human_PASdoublefilter.SAF -F SAF -o ../data/Test_FC_methods/Human_20 ../data/TestMM2/human*N.MM2.sort.bam  -s 2
