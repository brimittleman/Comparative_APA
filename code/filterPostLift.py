
def main(Orig,rev, outF):
    originalF=open(Orig,"r")
    revF=open(rev, "r")
    fout=open(outF,"w")
    okPAS={}
    for ln in revF:
        chrom, start, end, name, score, strand = ln.split()
        locID=chrom + ":" + start + ":" + "end"
        okPAS[name]=locID
    for ln in originalF:
        chrom, start, end, name, score, strand = ln.split()
        locID=chrom + ":" + start + ":" + "end"
        if name in okPAS.keys():
            if okPAS[name]==locID:
                fout.write(ln)
    fout.close()


if __name__ == "__main__":
    import sys
    inOrig =sys.argv[1]
    inRev= sys.argv[2]
    outFile= sys.argv[3]
    main(inOrig,inRev, outFile)
