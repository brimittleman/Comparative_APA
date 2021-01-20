def main(inFile, outFile):
    fout=open(outFile, "w")
    for ln in open(inFile, "r"):
        chr, start,end, name, score, strand =ln.split()
        if strand == "+":
            start_new= int(end)-100
            end_new=int(end)+100
        else:
            start_new= int(start)-100
            end_new=int(start)+100
        fout.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(chr, start_new, end_new, name, score ,strand))
    fout.close()


if __name__ == "__main__":
    import sys
    inFile =sys.argv[1]
    outFile= sys.argv[2]
    main(inFile, outFile)
