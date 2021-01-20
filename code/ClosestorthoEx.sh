#!/bin/bash

#SBATCH --job-name=ClosestorthoEx
#SBATCH --time=24:00:00
#SBATCH --partition=broadwl
#SBATCH --mem=32G
#SBATCH --mail-type=END
#SBATCH --output=ClosestorthoEx.out
#SBATCH --error=ClosestorthoEx.err

source ~/activate_anaconda.sh
conda activate comp_threeprime_env


bedtools closest -s -id -D "a" -a ../data/TestAnnoBiasOE/HumanIntronicGeneinOE.bed -b ../data/OrthoExonBed/human.noM.sort.merged.bed  > ../data/TestAnnoBiasOE/HumanUpstream.intronic.txt

bedtools closest -s -iu -D "a" -a ../data/TestAnnoBiasOE/HumanIntronicGeneinOE.bed -b ../data/OrthoExonBed/human.noM.sort.merged.bed > ../data/TestAnnoBiasOE/HumanDownstream.intronic.txt


bedtools closest -s -id -D "a" -a ../data/TestAnnoBiasOE/ChimpIntronicGeneinOE.bed -b ../data/OrthoExonBed/chimp.noM.sort.merged.bed  > ../data/TestAnnoBiasOE/ChimpUpstream.intronic.txt

bedtools closest -s -iu -D "a" -a ../data/TestAnnoBiasOE/ChimpIntronicGeneinOE.bed -b ../data/OrthoExonBed/chimp.noM.sort.merged.bed > ../data/TestAnnoBiasOE/ChimpDownstream.intronic.txt
