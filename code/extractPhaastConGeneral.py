def main(bedfile, fout):
  bw = pyBigWig.open("../data/PhastCon/hg38.phastCons100way.bw")
  region=open(bedfile,"r")
  fout=open(fout,"w")
  for i, line in enumerate(region):
      if i >1:
          #print(line)
          cols = line.strip().split()
          vals = bw.values(cols[0], int(cols[1]), int(cols[2]))
          valsNp=np.array(vals)
          meanVal=np.mean(valsNp)
          fout.write("%s\t%s\t%s\t%s\t%s\n"%(cols[0], cols[1], cols[2], cols[3], meanVal))
  fout.close()
  bw.close()


if __name__ == "__main__":
    import sys
    import numpy as np
    import pyBigWig
    inFile = sys.argv[1]
    outFile=sys.argv[2]
    main(inFile, outFile)
