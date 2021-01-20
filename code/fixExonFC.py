#/project2/gilad/briana/Comparative_APA/Human/data/RNAseq/sort/YG-BM-S10-18504H-Total_S10_R1_001-sort.bam

def main(infile, outfile):
    fout=open(outfile, "w")
    inf=open(infile,"r")
    for line, i in enumerate(inf):
        if line == 1:
            i_list=i.split()
            libraries=i_list[:6]
            for sample in i_list[6:]:
                full = sample.split("/")[9]
                print(full)
                samp= full.split("-")[3][:-1]
                samp_st="NA"+ samp
                libraries.append(samp_st)
            first_line= "\t".join(libraries)
            fout.write(first_line + '\n')
        else :
            fout.write(i)
    fout.close()


if __name__ == "__main__":
    from misc_helper import *
    import sys
    inFile =sys.argv[1]
    outFile= sys.argv[2]
    main(inFile, outFile)
