#!/bin/bash

#SBATCH --job-name=ALLPAS_sequenceDF
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=ALLPAS_sequenceDF.out
#SBATCH --error=ALLPAS_sequenceDF.err
#SBATCH --partition=broadwl
#SBATCH --mem=36G
#SBATCH --mail-type=END

source ~/activate_anaconda.sh
conda activate comp_threeprime_env

for i in AAAAAG AATACA AATAGA AATATA ACTAAA AGTAAA CATAAA GATAAA TATAAA AAAAAA  
do
echo $i
bedtools nuc -s  -pattern $i -C -fi /project2/gilad/kenneth/References/human/genome/hg38.fa -bed ../data/PAS_doubleFilter/PAS_doublefilter_either_HumanCoordHummanUsage.sort.bed > ../data/SignalSites_doublefilter/PAS_doublefilter_either_HumanCoordHummanUsage_${i}.txt
done 

for i in AAAAAG AATACA AATAGA AATATA ACTAAA AGTAAA CATAAA GATAAA TATAAA AAAAAA  
do
echo $i
bedtools nuc -s  -pattern $i  -C -fi /project2/gilad/briana/genome_anotation_data/Chimp_genome/panTro6.fa -bed ../data/PAS_doubleFilter/PAS_doublefilter_either_ChimpCoordChimpUsage.sort.bed > ../data/SignalSites_doublefilter/PAS_doublefilter_either_ChimpCoordChimpUsage_${i}.txt
done
