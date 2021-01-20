#this takes the meta file and finds genes with a dominant PAS above the cutoff percent, I can also give it a species to run the analysis on. The outfile will have the dominant PAS for the genes passing the criteria

def main(input, output, cutt, species):
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
    printdic={}
    # move through dic, find genes/pos dom
    outdic={}
    for key, value in PASdic.items():
        #print(key)
        #print(value)
        if len(value) >= 2:
            value.sort(reverse=True)
            #print(value)
            diffVal=value[0] - value[1]
            #print(diffVal)
            if diffVal >= float(cutt):
                geneUsage=key + ":" + str(value[0])
                outdic[geneUsage]= ""

    #write out
    #print(outdic)
    fin2=open(input, "r")
    fout=open(output, "w")
    for i, ln in enumerate(fin2):
        if i >0:
            if species == "Human":
                pas, disc, gene, loc, chr, start, end, Chimp, usage, Strand = ln.split()
                geneus=gene + ":" + usage
                #print(geneus)
                if gene not in printdic.keys():
                    if geneus in outdic.keys():
                        printdic[gene]=""
                        fout.write(ln)
            else:
                pas, disc, gene, loc, chr, start, end, usage, Human, Strand = ln.split()
                geneus=gene + ":" + usage
                if gene not in printdic.keys():
                    if geneus in outdic.keys():
                        printdic[gene]=""
                        fout.write(ln)
    fout.close()




if __name__ == "__main__":
    import sys
    inFile = sys.argv[1]
    outFile = sys.argv[2]
    cutoff = sys.argv[3]
    spec = sys.argv[4]
    main(inFile, outFile, cutoff, spec)
