#python file to name peaks


def main(inFile,spec, outFile):
    rows=0
    for ln in open(inFile, "r"):
        rows +=1
    rowname=list(range(0, rows))
    outF=open(outFile, "w")
    inF= open(inFile, "r")
    for i, ln in enumerate(inF):
        chrom, start, end, cov, strand, score = ln.split()
        #chromNoch=chrom[3:]
        name=spec + str(rowname[i])
        outF.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(chrom, start, end, name, cov, strand))
    outF.close()


if __name__ == "__main__":
    import numpy as np
    from misc_helper import *
    import sys
    inFile =sys.argv[1]
    species= sys.argv[2]
    outFile= sys.argv[3]
    main(inFile, species, outFile)
