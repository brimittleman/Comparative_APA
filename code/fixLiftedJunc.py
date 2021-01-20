def main(infile, outfile):
    inJunc= open(infile, "r")
    fout=open(outfile, "w")
    for ln in inJunc:
        line_list=ln.split()
        line_list[8]="255,0,0"
        line_list[-2]=line_list[-2][:-1]
        line_list[-1]=line_list[-1][:-1]
        #print(line_list)
        printLine="\t".join(line_list)
        fout.write("%s\n"%(printLine))
    fout.close()


if __name__ == "__main__":
  import sys
  inFile =sys.argv[1]
  outFile= sys.argv[2]
  main(inFile, outFile)
