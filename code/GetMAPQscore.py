def main(inFile, outFile):
    infile=open(inFile, "r")
    fout=open(outFile, "w")
    val={}
    for i in range(0, 256):
        chrval=str(i)
        val[chrval]=0
    print(val)
    for ln in infile:
        score1, score2, stat = ln.split()
        val[score2]+=1
    for key, val in val.items():
        fout.write("%s\t%d\n"%(key,val))
    fout.close()

if __name__ == "__main__":
    import sys
    inFile = sys.argv[1]
    outFile = sys.argv[2]
    main(inFile, outFile)
