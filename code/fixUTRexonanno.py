infile="/project2/gilad/briana/genome_anotation_data/Chimp_refseqAnno/pantro6_ncbiRefseq_Formatted_Allannotation.sort.bed"

outFile=open("/project2/gilad/briana/genome_anotation_data/Chimp_refseqAnno/pantro6_ncbiRefseq_Formatted_AllannotationUTRfix.sort.bed", "w")

#create a dictionary with the utr locations. only write out exons if the locations are not in the dictionary
utr_dic={}
for ln in open(infile,"r"):
    chrom, start, end, name, score, strand=ln.split()
    chrom_loc=chrom + ":" + start +":" + end
    loc= name.split(":")[0]
    gene=name.split(":")[1]
    if loc in ["utr3", "utr5"]:
        utr_dic[chrom_loc]=gene

for ln in open(infile, "r"):
  chrom, start, end, name, score, strand=ln.split()
  chrom_loc=chrom + ":" + start +":" + end
  loc= name.split(":")[0]
  gene=name.split(":")[1]
  if loc == "exon":
      if chrom_loc not in utr_dic.keys():
          outFile.write(ln)
  else:
    outFile.write(ln)
