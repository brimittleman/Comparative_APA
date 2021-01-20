def main(infile, outfile):
    okPAS={}
    okpasfile=open("../data/Peaks_5perc/Peaks_5perc_either_bothUsage_noUnchr.txt","r")
    for i, ln in enumerate(okpasfile):
        if i >0:
            pas=ln.split()[0]
            okPAS[pas]=""
        #print(okPAS)
    fin=open(infile,"r")
    fout=open(outfile,"w")
    for i, ln in enumerate(fin):
        if i ==0:
            fout.write(ln)
        else:
            #chr1:778450:778650:LOC100288069_+_utr3-Human-human31
            PAS= ln.split()[0].split("-")[-1]
            #print(PAS)
            if PAS in okPAS.keys():
                fout.write(ln)
    fout.close()


if __name__ == "__main__":
    import sys
    inFile = sys.argv[1]
    outFile = sys.argv[2]
    main(inFile, outFile)
