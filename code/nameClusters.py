#python file to name cluster

def main(inFile, outFile):
    rows=0
    for ln in open(inFile, "r"):
        rows +=1
    rowname=list(range(1, rows+1))
    outF=open(outFile, "w")
    inF= open(inFile, "r")
    for i, ln in enumerate(inF):
        chrom, start, end, name, score, strand = ln.split()
        score=rowname[i]
        outF.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(chrom, start, end, name, score, strand))
    outF.close()


if __name__ == "__main__":
    import numpy as np
    from misc_helper import *
    import sys
    inFile =sys.argv[1]
    outFile= sys.argv[2]
    main(inFile, outFile)
