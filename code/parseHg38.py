
txn=open("../../genome_anotation_data/hg38_refseq_anno/hg38_ncbiRefseq_Genes-2", "r")

fout=open("../../genome_anotation_data/hg38_refseq_anno/hg38_ncbiRefseq_GenesParsed.bed", "w")


geneDic={}
for i,ln in enumerate(txn):
    if i >0:
        tx, chrom, strand, start, end, geneN = ln.split()
        GeneChrom=geneN + "_" + chrom + "_" + strand
        if GeneChrom not in geneDic.keys():
            geneDic[GeneChrom]=[int(start),int(end)]
        else:
            geneDic[GeneChrom].append(int(start))
            geneDic[GeneChrom].append(int(end))

for key, values in geneDic.items():
    geneVec= key
    start= min(values)
    end=max(values)
    gene=geneVec.split("_")[0]
    chm= geneVec.split("_")[1]
    stran=geneVec.split("_")[2]
    score=0
    if stran in ["+", "-"]:
        fout.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(chm, start, end, gene, score, stran))
fout.close()
