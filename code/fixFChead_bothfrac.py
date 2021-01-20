def main(infile, outfile):
    fout=open(outfile, "w")
    inf=open(infile,"r")
    for line, i in enumerate(inf):
        if line == 1:
            i_list=i.split()
            libraries=i_list[:6]
            print(i_list)
            for sample in i_list[6:]:
                full = sample.split("/")[4]
                samp= full.split("-")[0].split("_")[2:]
                lim="_"
                samp_st=lim.join(samp)
                print(samp)
                libraries.append(samp_st)
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
