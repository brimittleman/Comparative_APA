
def main(inFile, outFile):
    fout=open(outFile,'w')
    fout.write("GeneID\tChr\tStart\tEnd\tStrand\n")
    for ln in open(inFile, "r"):
         chrom, start, end, name, score, strand = ln.split()
         ID = "peak_%s_%s_%s_%s"%(chrom,start, end, name)
         fout.write("%s\t%s\t%s\t%s\t%s\n"%(ID, chrom, start, end, strand))

    fout.close()


if __name__ == "__main__":
    import sys
    inFile =sys.argv[1]
    outFile= sys.argv[2]
    main(inFile, outFile)
