def main(inFile,outFile):
    #chrX	OrthoExon	exon	100632485	100632568	.	-	gene_id "ENSG00000000003"; gene_name "TSPAN6"; gene_version "15"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_id "ENST00000373020"; transcript_version "9"; transcript_name "TSPAN6-201"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; transcript_support_level "1"; exon_id "ENSE00000401072"; exon_number "6"; tag "basic";
    fin=open(inFile, "r")
    fout=open(outFile, "w")
    for i,ln in enumerate(fin):
        if i>=0:
            chrom=ln.split()[0]
            start=ln.split()[3]
            end=ln.split()[4]
            strand=ln.split()[6]
            #print(strand)
            #name is name and score is exon exon_number
            name=ln.split()[11][1:-2]
            #.split(";")[1][11:-1]
            #print(name)
            score=ln.split()[33][1:-2]
            if score == "NA":
                score=0
            #.split(";")[13][13:-1]
            #print(score)
            fout.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(chrom, start, end, name, score, strand))
    fout.close()

if __name__ == "__main__":
  import sys
  inSAF =sys.argv[1]
  outBed=sys.argv[2]
  main(inSAF, outBed)
