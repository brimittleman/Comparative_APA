#ortho UTR to fasta for ARE element

def main(genes, OutFile, ):
    seqin=open("../data/orthoUTR/HumanDistal3UTR.sort_withSeq.bed", "r")
    genedic={}
    fout=open(OutFile,"w")
    for ln in open(genes,"r"):
        g=ln.strip()
        genedic[g]=""
    for ln in seqin:
        gene=ln.split()[3]
        if gene in genedic.keys():
            fout.write(">%s\n"%gene)
            seq=ln.split()[15]
            fout.write("%s\n"%seq)
    fout.close()






if __name__ == "__main__":
    import sys
    geneList = sys.argv[1]
    Outfile =sys.argv[2]
    main(geneList, Outfile)
