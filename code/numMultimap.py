def main(inFile, outFile):
    infile=open(inFile, "r")
    fout=open(outFile, "w")
    nhval={}
    for i in range(1, 11):
        value= ":".join(["NH", "i", str(i)])
        nhval[value]=0
    for ln in infile:
        score1, score2, stat = ln.split()
        nhval[stat]+=1
    for key, val in nhval.items():
        fout.write("%s\t%d\n"%(key,val))
    fout.close()

if __name__ == "__main__":
    import sys
    inFile = sys.argv[1]
    outFile = sys.argv[2]
    main(inFile, outFile)
