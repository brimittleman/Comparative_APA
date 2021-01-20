#!/bin/bash

#SBATCH --job-name=JunctionLift
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=JunctionLift.out
#SBATCH --error=JunctionLift.err
#SBATCH --partition=broadwl
#SBATCH --mem=36G
#SBATCH --mail-type=END

source ~/activate_anaconda.sh
conda activate comp_threeprime_env

# oldFile map.chain newFile unMapped


echo "human 2 chimp"

for juncfile in `ls /project2/gilad/briana/Comparative_APA/Human/data/RNAseq/sort/*.junc`
do
    echo $juncfile
    liftOver $juncfile ../data/chainFiles/hg38ToPanTro6.over.chain $juncfile.2Chimp $juncfile.unlift -bedPlus=12
done


echo "human in chimp back to human"

for juncfile in `ls /project2/gilad/briana/Comparative_APA/Human/data/RNAseq/sort/*.junc.2Chimp`
do
    echo $juncfile
    liftOver $juncfile ../data/chainFiles/panTro6ToHg38.over.chain $juncfile.backHuman.junc $juncfile.unlift -bedPlus=12
done

#chimp, human, chimp
echo "chimp 2 human"

for juncfile in `ls /project2/gilad/briana/Comparative_APA/Chimp/data/RNAseq/sort/*.junc`
do
    echo $juncfile
    liftOver $juncfile ../data/chainFiles/panTro6ToHg38.over.chain $juncfile.2Human $juncfile.unlift -bedPlus=12
done

echo "chimp in human back to chimp"

for juncfile in `ls /project2/gilad/briana/Comparative_APA/Chimp/data/RNAseq/sort/*.junc.2Human`
do
    echo $juncfile
    liftOver $juncfile ../data/chainFiles/hg38ToPanTro6.over.chain  $juncfile.backChimp $juncfile.unlift -bedPlus=12
done
