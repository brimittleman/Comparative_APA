#!/bin/bash


#SBATCH --job-name=liftVCF
#SBATCH --output=liftVCF.out
#SBATCH --error=lliftVCF.err
#SBATCH --time=10:00:00
#SBATCH --partition=gilad
#SBATCH --nodelist=midway-l16b-31 
#SBATCH --mem=550G
#SBATCH --mail-type=END

module load picard 

#test chrom 4
java -jar $PICARD LiftoverVcf I=/project2/gilad/briana/li_genotypes/genotypesYRI.gen.proc.5MAF.chr4_test.vcf O=../data/QTLPASoverlap/Ch4_geno2pantro.vcf   CHAIN=../data/chainFiles/hg19ToPanTro6.over.chain.gz  REJECT=../data/QTLPASoverlap/Ch4_rejected_variants.vcf  R=/project2/gilad/briana/genome_anotation_data/Chimp_genome/panTro6.fa
