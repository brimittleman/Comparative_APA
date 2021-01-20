#> dummy1
#cagGTAAGT
#> dummy2
#gagGTAAGT
#> dummy3
#taaATAAGT

def main(infile, outfile):
  fout=open(outfile, "w")
  fin=open(infile, "r")
  for i, ln in enumerate(fin):
      if i >0:
          PAS=ln.split()[3]
          Seq=ln.split()[15]
          firstSeq=Seq[0:3].lower()
          endSeq=Seq[3:].upper()
          fullSeq=firstSeq + endSeq
          fout.write(">%s\n"%(PAS))
          fout.write("%s\n"%(fullSeq))
  fout.close()


if __name__ == "__main__":
  import numpy as np
  from misc_helper import *
  import sys
  inFile =sys.argv[1]
  outFile= sys.argv[2]
  main(inFile, outFile)
