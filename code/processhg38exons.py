#curate hg38 exons
#exon file: chr, start, end, strand,genename


exons={}
exonNames=open("../data/DiffSplice/hg38_ncbiRefseq_exonsNames", "r")
for ln in exonNames:
    exon, gene =ln.split()
    exons[exon]=gene


exonIn=open("../data/DiffSplice/hg38_ncbiRefseq_exons", "r")
exonOut=open("../data/DiffSplice/hg38_ncbiRefseq_exonsfixed", "w")
exonOut.write("chr\tstart\tend\tstrand\tgene_name\n")
for ln in exonIn:
    chr, start, end, name, score, strand= ln.split()
    exon=name.split("_")[0:2]
    exonfull="_".join(exon)
    gene=exons[exonfull]
    exonOut.write("%s\t%s\t%s\t%s\t%s\n"%(chr, start, end, strand, gene))

exonOut.close()
