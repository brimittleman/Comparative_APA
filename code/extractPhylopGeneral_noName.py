def main(bedfile, fout):
  bw = pyBigWig.open("../data/PhyloP/hg38.phyloP100way.bw")
  region=open(bedfile,"r")
  fout=open(fout,"w")
  for i, line in enumerate(region):
      if i >1:
          #print(line)
          cols = line.strip().split()
          if cols[0] in ["chr1", "chr2", "chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"]:
              vals = bw.values(cols[0], int(cols[1]), int(cols[2]))
              valsNp=np.array(vals)
              meanVal=np.mean(valsNp)
              fout.write("%s\t%s\t%s\t%s\n"%(cols[0], cols[1], cols[2], meanVal))
  fout.close()
  bw.close()


if __name__ == "__main__":
    import sys
    import numpy as np
    import pyBigWig
    inFile = sys.argv[1]
    outFile=sys.argv[2]
    main(inFile, outFile)
