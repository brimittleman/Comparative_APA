

def main(input, output, species):
    PASdic={}
    fin=open(input, "r")
    #human27 Both LOC105378580 intron chr1 788172 788372 0.331666666666667 0.13 +
     #create dic
    for i, ln in enumerate(fin):
        if i >0:
            if species == "Human":
                pas, disc, gene, loc, chr, start, end, Chimp, usage, Strand = ln.split()
                if float(usage) > 0:
                    if gene not in PASdic.keys():
                            PASdic[gene]=[float(usage)]
                    else:
                        PASdic[gene].append(float(usage))
            else:
                pas, disc, gene, loc, chr, start, end, usage, Human, Strand = ln.split()
                if float(usage) > 0:
                    if gene not in PASdic.keys():
                        PASdic[gene]=[float(usage)]
                    else:
                        PASdic[gene].append(float(usage))
    fin.close()
    #print(PASdic)
    # move through dic, find genes/pos dom
    outdic={}
    for key, value in PASdic.items():
        #print(key)
        #print(value)
        if len(value) >= 2:
            value.sort(reverse=True)
            #print(value)
            diffVal=value[0] - value[1]
            geneUsage=key + ":" + str(value[0])
            outdic[geneUsage]= diffVal

    #write out
    #print(outdic)
    printdic={}
    fin2=open(input, "r")
    fout=open(output, "w")
    for i, ln in enumerate(fin2):
        if i >0:
            if species == "Human":
                pas, disc, gene, loc, chr, start, end, Chimp, usage, Strand = ln.split()
                geneus=gene + ":" + usage
                if gene not in printdic.keys():
                    if geneus in outdic.keys():
                        diff=outdic[geneus]
                        printdic[gene]=""
                        fout.write("%s\t%s\t%s\n"%(pas, gene, diff))
            else:
                pas, disc, gene, loc, chr, start, end, usage, Human, Strand = ln.split()
                geneus=gene + ":" + usage
                if gene not in printdic.keys():
                    if geneus in outdic.keys():
                        diff=outdic[geneus]
                        printdic[gene]=""
                        fout.write("%s\t%s\t%s\n"%(pas, gene, diff))
    fout.close()




if __name__ == "__main__":
    import sys
    inFile = sys.argv[1]
    outFile = sys.argv[2]
    spec = sys.argv[3]
    main(inFile, outFile, spec)
