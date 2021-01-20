def main(infile, outfile):
    inJunc= open(infile, "r")
    fout=open(outfile, "w")

    chr_dic={}
    for i in range(1,23):
        toadd= "chr" + str(i)
        chr_dic[toadd]=""
    chr_dic["chr2A"]=""
    chr_dic["chr2B"]=""
    chr_dic["X"]=""
    chr_dic["Y"]=""
    for ln in inJunc:
        chrom=ln.split()[0]
        if chrom in chr_dic:
            line_list=ln.split()
            line_list[8]="255,0,0"
            line_list[-2]=line_list[-2][:-1]
            line_list[-1]=line_list[-1][:-1]
            printLine="\t".join(line_list)
            fout.write("%s\n"%(printLine))
    fout.close()


if __name__ == "__main__":
  import sys
  inFile =sys.argv[1]
  outFile= inFile + ".fixed.junc"
  main(inFile, outFile)
