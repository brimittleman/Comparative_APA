def main(inFile, outFile):
    fin=open(inFile, "r")
    fout=open(outFile, "w")
    for ln in fin:
        ln_listpos=ln.split()
        fout.write(ln)
        ln_listpos[-2]= "-"
        newLn="\t".join(ln_listpos)
        fout.write("%s\n"%(newLn))
    fout.close()


if __name__ == "__main__":
    import sys
    inFile = sys.argv[1]
    outFile = sys.argv[2]
    main(inFile, outFile)
