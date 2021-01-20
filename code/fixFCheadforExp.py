def main(infile, outfile):
    fout=open(outfile, "w")
    inf=open(infile,"r")
    #../data/TestMM2/human_combined_18523_N.MM2.sort.bam
    for line, i in enumerate(inf):
        if line == 1:
            i_list=i.split()
            libraries=i_list[:6]
            for sample in i_list[6:]:
                full = sample.split("/")[3]
                samp= full.split(".")[0]
                #amp_st=lim.join(samp)
                libraries.append(samp)
            first_line= "\t".join(libraries)
            fout.write(first_line + '\n')
        else :
            fout.write(i)
    fout.close()


if __name__ == "__main__":
    import numpy as np
    from misc_helper import *
    import sys
    inFile =sys.argv[1]
    outFile= sys.argv[2]
    main(inFile, outFile)
