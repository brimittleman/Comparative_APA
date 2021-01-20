TXN2Gene_file="/project2/gilad/briana/genome_anotation_data/hg38_refseq_anno/hg38_ncbiRefseq_txn2genename"

gene_dic={}

for ln in open(TXN2Gene_file,"r"):
   txn=ln.split()[1]
   gene=ln.split()[12]
   gene_dic[txn]=gene

outF=open("/project2/gilad/briana/genome_anotation_data/hg38_refseq_anno/hg38_ncbiRefseq_Formatted_Allannotation.sort.bed","w")

inFile="/project2/gilad/briana/genome_anotation_data/hg38_refseq_anno/hg38_ncbiRefseq_Allannotation.sort.bed"  


for ln in open(inFile, "r"):
   chrom, start, end, name, score, strand = ln.split()
   chrom_fix=chrom[3:]
   txn=name.split("_")[:2]
   #print(txn)
   txn2="_".join(txn)
   gene=gene_dic[txn2]
   anno=name.split("_")[2]
   ids =anno + ":"+ gene
   outF.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(chrom, start, end, ids, score, strand))

outF.close()
